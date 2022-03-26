import os
import csv

import networkx as nx

from scripts.utils import make_graph_name


RANKS = [
    "kingdom",
    "supergroup",
    "division",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


MAX_EDGES = 50000


def filter_graph(graph: nx.Graph, n_edges: int) -> nx.Graph:
    new_graph = nx.Graph()
    pos_edges = 0
    neg_edges = 0
    for _, _, attrs in graph.edges(data=True):
        if attrs["sign"] == "+":
            pos_edges += 1
        if attrs["sign"] == "-":
            neg_edges += 1
    try:
        pos_to_keep = int(n_edges * pos_edges / (pos_edges + neg_edges))
        neg_to_keep = int(n_edges * neg_edges / (pos_edges + neg_edges))
    except ZeroDivisionError as e:
        print(pos_edges, neg_edges)
        print(graph)
        pos_to_keep = n_edges
        neg_to_keep = 0
    print(f"Keeping {pos_to_keep} positive edges and {neg_to_keep} negative edges")
    for head, tail, attrs in sorted(
        graph.edges(data=True), key=lambda x: x[2]["weight"], reverse=True
    )[:pos_to_keep]:
        new_graph.add_edge(head, tail, **attrs)
    for head, tail, attrs in sorted(
        graph.edges(data=True), key=lambda x: -x[2]["weight"], reverse=True
    )[:neg_to_keep]:
        new_graph.add_edge(head, tail, **attrs)
    return new_graph


def parse_edgelist(fname, delimiter):
    if delimiter is not None:
        return nx.read_edgelist(fname, data=(("weight", float),), delimiter="\t")
    else:
        return nx.read_edgelist(fname, data=(("weight", float),))


def read_tsv(fname):
    graph = nx.Graph()
    with open(fname) as handle:
        fl = next(handle)
        if fl.strip() == '"from"  "to"':
            for line in handle:
                _, left, right = line.strip().split("\t")
                graph.add_edge(
                    left.strip('"').lstrip("X"), right.strip('"').lstrip("X")
                )
        else:
            data = []
            row_names = []
            for line in handle:
                row_name, *d = line.strip().split("\t")
                row_names.append(row_name)
                data.append(d)
            for left, row in zip(row_names, data):
                for right, value in zip(row_names, row):
                    if value == "1":
                        graph.add_edge(left.strip('"'), right.strip('"'))
    return graph


def read_fastspar(fname):
    graph = nx.Graph()
    config = snakemake.config["fastspar_configs"][os.path.split(fname)[1]]
    val_type = config.get("network_base", "correlations")
    with open(os.path.join(fname, f"{val_type}.tsv")) as handle:
        reader = csv.reader(handle, delimiter="\t")
        row_names = next(reader)[1:]
        for left, *row in reader:
            for right, value in zip(row_names, row):
                if abs(float(value)) > config["threshold"] and left != right:
                    graph.add_edge(left, right, weight=float(value))
    return graph


def parse_conet(fpath):
    graph = nx.Graph()
    with open(fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for line in reader:
            # for now just getting weight as spearman correlation
            row = dict(zip(header, line))
            graph.add_edge(
                row["head"], row["tail"], weight=float(row["spearman_weight"])
            )
    return graph


for filepath in snakemake.input["networks"]:
    if os.path.split(filepath)[0].endswith("conet_results"):
        graph = parse_conet(filepath)
    elif filepath.endswith(".edgelist"):
        if "flashweave" in filepath:
            delimiter = "\t"
        else:
            delimiter = None
        graph = parse_edgelist(filepath, delimiter)
    elif filepath.endswith(".tsv"):
        graph = read_tsv(filepath)
    elif os.path.isdir(filepath):
        graph = read_fastspar(filepath)
    else:
        raise ValueError("Unrecognized network format!")
    for head, tail, attrs in graph.edges(data=True):
        graph[head][tail]["sign"] = "+" if attrs["weight"] >= 0 else "-"
    graph = filter_graph(graph, MAX_EDGES)
    nx.write_edgelist(
        graph, make_graph_name(filepath), data=["weight", "sign"], delimiter="\t"
    )