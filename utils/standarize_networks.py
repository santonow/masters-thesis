import os
import csv

import networkx as nx

from occurence_table import OTUTable
from .utils import make_graph_name

table = OTUTable.from_tsv(snakemake.input["base"], snakemake.input["meta"])


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


for filepath in snakemake.input['networks']:
    if filepath.endswith('.edgelist'):
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
    for _id in table.table.ids(axis='observation'):
        graph.add_node(_id)
    nx.write_edgelist(graph, make_graph_name(filepath))
