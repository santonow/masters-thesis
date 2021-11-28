from collections import defaultdict
import os

import networkx as nx
import numpy as np


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


class Indexer:

    def __init__(self):
        self.name_to_id = dict()

    def __call__(self, s):
        if s not in self.name_to_id:
            self.name_to_id[s] = len(self.name_to_id)
        return self.name_to_id[s]


def prepare_files(graph_path, taxonomy_file, max_taxa=100, tax_level=None):
    if max_taxa is None:
        max_taxa = 100
    if tax_level is None:
        lineage_min_length = -1
    else:
        lineage_min_length = RANKS.index(tax_level)
    net = nx.read_edgelist(
        graph_path,
        delimiter='\t',
        data=(('weight', float),)
    )

    max_edges = max_taxa * 3

    with open(taxonomy_file) as handle:
        import csv
        taxonomy = dict()
        reader = csv.reader(handle, delimiter='\t')
        for line in reader:
            taxonomy[line[0]] = tuple(line[1:])

    edges_per_node = defaultdict(lambda: defaultdict(list))

    for head, tail, attrs in net.edges(data=True):
        if head in taxonomy and tail in taxonomy:
            if taxonomy[head] != taxonomy[tail]:
                if len(taxonomy[head]) > lineage_min_length and len(taxonomy[tail]) > lineage_min_length:
                    edges_per_node[taxonomy[head]][taxonomy[tail]].append(attrs['weight'])

    filtered_graph = nx.Graph()
    for head, neighbors in edges_per_node.items():
        for tail, weights in neighbors.items():
            if len(weights) > 1:
                filtered_graph.add_edge(head, tail, weight=len(weights))

    indexer = Indexer()
    hash_to_tax = dict()

    for head, tail in list(filtered_graph.edges()):
        if head == tail and filtered_graph.has_edge(head, tail):
            filtered_graph.remove_edge(head, tail)
    for node in list(filtered_graph.nodes()):
        if node == ('Eukaryota',):
            filtered_graph.remove_node(node)
    for node, degree in list(filtered_graph.degree()):  # noqa
        if degree == 0:
            filtered_graph.remove_node(node)

    taxa = set()
    new_graph = nx.Graph()
    for head, tail, attrs in sorted(
        filtered_graph.edges(data=True), reverse=True, key=lambda x: abs(x[2]['weight'])
    ):
        head_hash = indexer(head)
        hash_to_tax[head_hash] = head
        tail_hash = indexer(tail)
        hash_to_tax[tail_hash] = tail
        if new_graph.number_of_edges() >= max_edges:
            break
        if len(taxa) < max_taxa:
            new_graph.add_edge(head_hash, tail_hash, **attrs)
            taxa.add(head_hash)
            taxa.add(tail_hash)
        elif head_hash in taxa and tail_hash in taxa:
            new_graph.add_edge(head_hash, tail_hash, **attrs)

    fname = os.path.split(graph_path)[1].split('.')[0]
    if tax_level is not None:
        new_dir = os.path.join(os.path.split(graph_path)[0], tax_level + '_' + fname)
    else:
        new_dir = os.path.join(os.path.split(graph_path)[0], fname)
    os.makedirs(new_dir, exist_ok=True)
    with open(os.path.join(new_dir, 'reduced_taxonomy.tsv'), 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['otu_id'] + RANKS)
        for h, lineage in hash_to_tax.items():
            if h in new_graph:
                new_lineage = list(lineage) + [lineage[-1] for _ in range(8 - len(lineage))]
                writer.writerow([h] + list(new_lineage))

    with open(os.path.join(new_dir, 'reduced_graph.edgelist'), 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for head, tail, attrs in new_graph.edges(data=True):
            writer.writerow([head, tail, attrs['weight']])




