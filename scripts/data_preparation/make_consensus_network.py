from collections import defaultdict
import os
from collections import Counter

import networkx as nx
import numpy as np


def group_by_program(filenames):
    result = defaultdict(list)
    for filename in filenames:
        dirname, _ = os.path.split(filename)
        result[dirname].append(
            nx.read_edgelist(
                filename, delimiter="\t", data=(("weight", float), ("sign", str))
            )
        )
    return result


consensus_network = nx.Graph()
edges_with_sources = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for dirname, networks in group_by_program(snakemake.input).items():
    for i, network in enumerate(networks):
        for head, tail, attrs in network.edges(data=True):
            edges_with_sources[tuple(sorted((head, tail)))][dirname][i].append(
                (attrs["weight"], attrs["sign"])
            )

for edge, edges_per_method in edges_with_sources.items():
    mean_weights = []
    signs = []
    for dirname, method in edges_per_method.items():
        method_weight = []
        for _, weights in method.items():
            method_weight.append(np.mean([x[0] for x in weights]))
            signs.extend(x[1] for x in weights)
        mean_weights.append(np.mean(method_weight))
    signs = Counter(signs)
    if signs["+"] == signs["-"]:
        sign = "?"
    else:
        sign = signs.most_common()[0][0]
    if len(edges_per_method) >= snakemake.config["consensus_network"]["min_programs"]:
        consensus_network.add_edge(*edge, weight=np.mean(mean_weights), sign=sign)

nx.write_edgelist(
    consensus_network, snakemake.output[0], data=["weight", "sign"], delimiter="\t"
)