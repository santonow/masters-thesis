import os
from collections import Counter
from functools import partial
from typing import Dict, Callable, Union, Any
import csv
import string

import networkx as nx
import numpy as np
import openpyxl
from joblib import Parallel, delayed


class OpenPyXLHandler:

    def __init__(self, fpath, sheet_name=None):
        self.unnamed_sheet_count = 0
        self.row_counts = Counter()
        self.fpath = fpath
        self.wb = openpyxl.Workbook()
        self.ws = self.wb.active
        if sheet_name is not None:
            self.ws.title = sheet_name
        else:
            self.ws.title = f'Sheet {self.unnamed_sheet_count}'
            self.unnamed_sheet_count += 1
        self.row_counts[self.ws.title] = 1

    def writerow(self, row):
        current_indice = self.row_counts[self.ws.title]
        for char, elem in zip(string.ascii_uppercase, row):
            self.ws[f'{char}{current_indice}'] = elem
        self.row_counts[self.ws.title] += 1

    def save(self):
        self.wb.save(self.fpath)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()


def get_n_signs(graph: nx.Graph, sign: str) -> int:
    n = 0
    for _, _, attrs in graph.edges(data=True):
        if attrs['sign'] == sign:
            n += 1
    return n


def return_if_except(exception, value: Any):
    def decorator(fun):
        def wrapper(*args, **kwargs):
            try:
                return fun(*args, **kwargs)
            except exception:
                return value
        return wrapper
    return decorator


@return_if_except(nx.NetworkXError, np.inf)
def get_radius(graph: nx.Graph) -> Union[int, float]:
    return nx.radius(graph)


@return_if_except(nx.NetworkXError, np.inf)
def get_diameter(graph: nx.Graph) -> Union[int, float]:
    return nx.diameter(graph)


METRIC_TO_FUN: Dict[str, Callable[[nx.Graph], Union[float, int]]] = {
    'Number of nodes': lambda graph: graph.number_of_nodes(),
    'Number of edges': lambda graph: graph.number_of_edges(),
    'Number of positive interactions': partial(get_n_signs, sign='+'),
    'Number of negative interactions': partial(get_n_signs, sign='-'),
    'Number of connected components': lambda graph: nx.number_connected_components(graph),
    'Mean network centrality': lambda graph: np.mean(list(nx.degree_centrality(graph).values())),
    'Mean network betweeness': lambda graph: np.mean(list(nx.betweenness_centrality(graph).values())),
    'Network radius': lambda graph: nx.radius(graph),
    'Network diameter': lambda graph: nx.diameter(graph),
    'Average clustering coefficient': lambda graph: nx.average_clustering(graph),
    'Average path length': lambda graph: nx.average_shortest_path_length(graph),
}
HEADER = list(METRIC_TO_FUN.keys())

METHOD_CAPITALIZATION = {
    'fastspar_results': 'FastSpar (SparCC)',
    'flashweave_results': 'FlashWeave',
    'phyloseq_results': 'phyloseq',
    'conet_results': 'Custom CoNet'
}


def prepare_stats(graph: nx.Graph):
    return {
        name: fun(graph) for name, fun in METRIC_TO_FUN.items()
    }


def read_graph(fpath: str) -> nx.Graph:
    graph = nx.Graph()
    with open(fpath) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for node1, node2, *attrs in reader:
            graph.add_edge(node1, node2, sign=attrs[1])
    return graph


def process(fpath: str):
    graph = read_graph(fpath)
    return fpath, prepare_stats(graph)


if __name__ == '__main__':
    proc_delayed = delayed(process)
    executor = Parallel(n_jobs=snakemake.threads, verbose=5)
    with OpenPyXLHandler(snakemake.output[0], 'Network statistics') as handle:
        handle.writerow(['network'] + HEADER)
        tasks = (proc_delayed(fpath) for fpath in snakemake.input)
        for fpath, stats in executor(tasks):
            if 'consensus_network' in fpath:
                network_name = 'Consensus network'
            else:
                head, tail = os.path.split(fpath)
                rank = tail.split('_')[0]
                method_name = METHOD_CAPITALIZATION[os.path.split(head)[1]]
                network_name = method_name + f' ({rank}) network'
            handle.writerow([network_name] + [stats[statistic] for statistic in HEADER])
