import csv
import json
from math import isnan
from collections import defaultdict, Counter
from typing import Dict, Tuple, Set

import pandas as pd
import numpy as np
import networkx as nx
import cartoGRAPHs


COLORS = [
    '#007399',
    '#008060',
    '#00b300',
    '#999900',
    '#b35900',
    '#ff0055',
    '#8a428a',
    '#1f42ad',
    '#7de8e8',
    '#bf80ff'
]


sources: Dict[str, str] = {
    k: v for k, v in snakemake.input.items()
    if k not in {'taxonomy', 'trophic_groups'}
}


def read_graph(
    sources: Dict[str, str], otu_to_lineage: Dict[str, Tuple[str, ...]], taxa_to_trophic_group: Dict[str, Set[str]]
) -> nx.Graph:
    graph = nx.Graph()
    for method, fpath in sources.items():
        with open(fpath) as handle:
            reader = csv.reader(handle, delimiter='\t')
            for head, tail, weight, sign in reader:
                if not graph.has_edge(head, tail):
                    graph.add_edge(head, tail)
                    graph[head][tail]['sources'] = []
                    graph[head][tail]['weights'] = []
                    graph[head][tail]['signs'] = []
                graph[head][tail]['sources'].append(method)
                graph[head][tail]['weights'].append(float(weight))
                graph[head][tail]['signs'].append(sign)

    for head, tail, attrs in graph.edges(data=True):
        weights, signs = [], []
        for source, weight, sign in zip(attrs['sources'], attrs['weights'], attrs['signs']):
            if source == 'consensus':
                graph[head][tail]['weight'] = weight
                graph[head][tail]['sign'] = sign
                break
            else:
                weights.append(weight)
                signs.append(sign)
        else:
            graph[head][tail]['weight'] = np.mean(weights)
            if signs.count('+') > signs.count('-'):
                graph[head][tail]['sign'] = '+'
            elif signs.count('+') < signs.count('-'):
                graph[head][tail]['sign'] = '-'
            else:
                graph[head][tail]['sign'] = '?'

    for node in graph.nodes():
        graph.nodes[node]['lineage'] = otu_to_lineage.get(node, tuple())
        for trophic_group, taxa in taxa_to_trophic_group.items():
            if any(elem in taxa for elem in graph.nodes[node]['lineage']):
                graph.nodes[node]['trophic_group'] = trophic_group
                break
        else:
            graph.nodes[node]['trophic_group'] = ''

    return graph


def collect_statistics(graph, all_trophic_groups: Set[str]):
    stats = dict()
    stats['degree centrality'] = nx.degree_centrality(graph)
    stats['betweenness centrality'] = nx.betweenness_centrality(graph)
    stats['eigen centrality'] = nx.eigenvector_centrality(graph)

    signs = defaultdict(Counter)
    weights = defaultdict(list)

    for head, tail, attrs in graph.edges(data=True):
        for node in [head, tail]:
            signs[node][attrs['sign']] += 1
            weights[node].append(abs(attrs['weight']))
    stats['pos'] = {
        node: c['+'] / nx.degree(graph, node) for node, c in signs.items()
    }
    stats['neg'] = {
        node: c['-'] / nx.degree(graph, node) for node, c in signs.items()
    }
    stats['weights'] = {
        node: np.mean(ws) for node, ws in weights.items()
    }
    for trophic_group in all_trophic_groups:
        stats[trophic_group] = Counter()
        for node in graph.nodes():
            if trophic_group == graph.nodes[node]['trophic_group']:
                stats[trophic_group][node] = 1
            else:
                stats[trophic_group][node] = 0

    print('Done.')
    return stats


def load_to_pandas(stats):
    node_stats = defaultdict(list)
    keys = sorted(stats)
    for node in stats['eigen centrality']:
        for key in keys:
            node_stats[node].append(stats[key][node])
    stats_df = pd.DataFrame.from_dict(node_stats, orient='index', columns=keys)
    return stats_df


def determine_sizes(stats, field='degree centrality'):
    ranges = dict()

    for i, elem in enumerate(np.array_split(sorted(list(stats[field].values())), 9), 1):
        print(len(elem))
        ranges[(elem[0], elem[-1])] = i / 2 + 1

    return ranges


def assign_range(node, sizes):
    node_centrality = stats['degree centrality'][node]
    for (left, right), i in sizes.items():
        if left <= node_centrality <= right:
            return i
    return 1


if __name__ == '__main__':
    trophic_groups = pd.read_excel(snakemake.input['trophic_groups'])
    taxa_to_trophic_group = dict()
    for row in trophic_groups.to_dict(orient='records'):
        if isinstance(row['Trophic group'], float) and isnan(row['Trophic group']):
            continue
        else:
            for i in reversed(range(1, 6)):
                if not (isinstance(row[f'Taxa_{i}'], float) and isnan(row[f'Taxa_{i}'])):
                    taxa_to_trophic_group[row[f'Taxa_{i}'].strip()] = row['Trophic group'].strip().replace(
                        'parasitic', 'parasite'
                    )
    group_to_taxa = defaultdict(set)
    for taxa, group in taxa_to_trophic_group.items():
        group_to_taxa[group].add(taxa)

    otu_to_lineage = dict()
    with open(snakemake.input['taxonomy']) as handle:
        reader = csv.reader(handle, delimiter='\t')
        print(next(reader))
        for line in reader:
            otu_to_lineage[line[0]] = tuple(line[1:])

    graph = read_graph(sources, otu_to_lineage, group_to_taxa)
    stats = collect_statistics(graph, set(group_to_taxa.keys()))

    layout = cartoGRAPHs.generate_layout(
        graph,
        dim=2,
        layoutmethod='precalculated',
        dimred_method='umap',
        Matrix=load_to_pandas(stats),
    )

    sizes = determine_sizes(stats)

    for node in list(graph.nodes()):
        if node in layout:
            graph.nodes[node]['x'], graph.nodes[node]['y'] = layout[node]

    higher_taxa = {x[1] if len(x) > 1 else x[0] for x in otu_to_lineage.values()}
    taxon_to_color = dict(zip(sorted(higher_taxa), COLORS))

    output_dict = {
        'nodes': [],
        'edges': []
    }

    for i, node in enumerate(graph.nodes()):
        if node != 'leftover_vector':
            lineage = otu_to_lineage[str(node)]
            if len(lineage) > 1:
                color = taxon_to_color[lineage[1]]
            else:
                color = taxon_to_color[lineage[0]]
            d = {
                'key': str(node),
                'attributes': {
                    'x': graph.nodes[node]['x'],
                    'y': graph.nodes[node]['y'],
                    'label': f'{node} ({", ".join(otu_to_lineage[str(node)][-3:])})',
                    'lineage': otu_to_lineage[str(node)],
                    'color': color,
                    'size': assign_range(str(node), sizes),
                }
            }
            output_dict['nodes'].append(d)

    for i, (head, tail, attrs) in enumerate(graph.edges(data=True)):

        if head != 'leftover_vector' and tail != 'leftover_vector':
            color = (128, 128, 128, 0.2)
            if attrs['sign'] == '+':
                color = (102, 102, 255, 0.2)
            elif attrs['sign'] == '-':
                color = (255, 153, 153, 0.2)
            d = {
                'key': str(i),
                'source': head,
                'target': tail,
                'attributes': {
                    'size': 0.5,
                    # 'color': color,
                    'sources': attrs['sources'],
                    'weight': attrs['weight'],
                    'sign': attrs['sign']
                }
            }
            output_dict['edges'].append(d)

    with open(snakemake.output[0], 'w') as handle:
        json.dump(output_dict, handle)
