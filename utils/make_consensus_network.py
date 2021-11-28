from collections import defaultdict
import os

import networkx as nx
import numpy as np
from plot.prepare_vis import prepare_files


def group_by_program(filenames):
    result = defaultdict(list)
    for filename in filenames:
        dirname, _ = os.path.split(filename)
        result[dirname].append(
            nx.read_edgelist(filename, delimiter='\t', data=(('weight', float),))
        )
    return result


consensus_network = nx.Graph()
edges_with_sources = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for dirname, networks in group_by_program(snakemake.input[:-1]).items():
    for i, network in enumerate(networks):
        for head, tail, attrs in network.edges(data=True):
            edges_with_sources[(head, tail)][dirname][i].append(attrs['weight'])

for edge, edges_per_method in edges_with_sources.items():
    mean_weights = []
    for dirname, method in edges_per_method.items():
        method_weight = []
        for _, weights in method.items():
            method_weight.append(np.mean(weights))
        mean_weights.append(np.mean(method_weight))
    if len(edges_per_method) >= snakemake.config['consensus_network']['min_programs']:
        consensus_network.add_edge(*edge, weight=np.mean(mean_weights))

nx.write_edgelist(
    consensus_network, snakemake.output[0], data=['weight'], delimiter='\t'
)

prepare_files(
    snakemake.output[0], snakemake.input[-1], snakemake.config['visualization']['max_nodes']
)
