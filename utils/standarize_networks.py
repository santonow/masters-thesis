import os
import csv

import networkx as nx

from utils import make_graph_name
from plot.prepare_vis import prepare_files


RANKS = [
    'kingdom',
    'supergroup',
    'division',
    'class',
    'order',
    'family',
    'genus',
    'species'
]


MAX_EDGES = 50000


def filter_graph(graph: nx.Graph, n_edges: int) -> nx.Graph:
    new_graph = nx.Graph()
    for head, tail, attrs in sorted(
        graph.edges(data=True), key=lambda x: abs(x[2]['weight']), reverse=True
    )[:n_edges]:
        new_graph.add_edge(head, tail, **attrs)
    return new_graph


def parse_edgelist(fname, delimiter):
    if delimiter is not None:
        return nx.read_edgelist(fname, data=(('weight', float),), delimiter='\t')
    else:
        return nx.read_edgelist(fname, data=(('weight', float),))


def read_tsv(fname):
    graph = nx.Graph()
    with open(fname) as handle:
        fl = next(handle)
        if fl.strip() == '"from"  "to"':
            for line in handle:
                _, left, right = line.strip().split('\t')
                graph.add_edge(
                    left.strip('"').lstrip('X'), right.strip('"').lstrip('X')
                )
        else:
            data = []
            row_names = []
            for line in handle:
                row_name, *d = line.strip().split('\t')
                row_names.append(row_name)
                data.append(d)
            for left, row in zip(row_names, data):
                for right, value in zip(row_names, row):
                    if value == '1':
                        graph.add_edge(left.strip('"'), right.strip('"'))
    return graph


def read_fastspar(fname):
    graph = nx.Graph()
    config = snakemake.config['fastspar_configs'][os.path.split(fname)[1]]
    val_type = config.get('network_base', 'correlations')
    with open(os.path.join(fname, f'{val_type}.tsv')) as handle:
        reader = csv.reader(handle, delimiter='\t')
        row_names = next(reader)[1:]
        for left, *row in reader:
            for right, value in zip(row_names, row):
                if float(value) > config['threshold'] and left != right:
                    graph.add_edge(left, right, weight=float(value))
    return graph


def parse_conet(fpath):
    graph = nx.Graph()
    with open(fpath) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        for line in reader:
            # for now just getting weight as spearman correlation
            row = dict(zip(header, line))
            graph.add_edge(row['head'], row['tail'], weight=float(row['spearman_weight']))
    return graph


for filepath in snakemake.input['networks']:
    if os.path.split(filepath)[0].endswith('conet_results'):
        graph = parse_conet(filepath)
    elif filepath.endswith('.edgelist'):
        if 'flashweave' in filepath:
            delimiter = '\t'
        else:
            delimiter = None
        graph = parse_edgelist(filepath, delimiter)
    elif filepath.endswith('.tsv'):
        graph = read_tsv(filepath)
    elif os.path.isdir(filepath):
        graph = read_fastspar(filepath)
    else:
        raise ValueError('Unrecognized network format!')
    graph = filter_graph(graph, MAX_EDGES)
    for head, tail, attrs in graph.edges(data=True):
        graph[head][tail]['sign'] = '+' if attrs['weight'] >= 0 else '-'
    nx.write_edgelist(
        graph, make_graph_name(filepath), data=['weight', 'sign'], delimiter='\t'
    )
    for rank in RANKS:
        prepare_files(
            make_graph_name(filepath),
            snakemake.input['tax_table'],
            snakemake.config.get('visualization', dict()).get('max_nodes'),
            tax_level=rank
        )
    prepare_files(
        make_graph_name(filepath),
        snakemake.input['tax_table'],
        snakemake.config.get('visualization', dict()).get('max_nodes'),
        tax_level=None
    )
