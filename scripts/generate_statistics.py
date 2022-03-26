import json
import os
from collections import defaultdict
from functools import partial
from itertools import combinations, chain
from typing import Dict, Callable, Union, Any, List, Optional, Tuple
import csv
import string
from math import isnan

import networkx as nx
import numpy as np
import openpyxl
import pandas as pd
from joblib import Parallel, delayed


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


class OpenPyXLWriter:
    def __init__(self, fname: str, sheet_name=None):
        self.fname = fname
        self.wb = openpyxl.Workbook()
        self.ws = self.wb.active
        if sheet_name is not None:
            self.ws.title = sheet_name
        self.row_count = {self.ws.title: 1}

    def _add_sheet(self, sheet_name: str):
        self.wb.create_sheet(title=sheet_name)
        self.row_count[sheet_name] = 1

    def _switch_sheet(self, sheet_name: str):
        if sheet_name not in self.row_count:
            self._add_sheet(sheet_name)
        self.ws = self.wb[sheet_name]

    def _update_row(self):
        self.row_count[self.ws.title] += 1

    @property
    def _current_pos(self):
        return self.row_count[self.ws.title]

    def _pos_iter(self):
        for letter in string.ascii_uppercase:
            yield f"{letter}{self._current_pos}"

    def writerow(self, row: List[Any], sheet_name: Optional[str] = None):
        if sheet_name is not None and self.ws.title != sheet_name:
            self._switch_sheet(sheet_name)
        for elem, pos in zip(row, self._pos_iter()):
            try:
                if elem is None:
                    self.ws[pos] = str(elem)
                else:
                    self.ws[pos] = elem
            except ValueError as e:
                self.ws[pos] = str(elem)
        self._update_row()

    def save(self):
        for sheet_name, count in self.row_count.items():
            if count == 1:
                del self.wb[sheet_name]
        self.wb.save(self.fname)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()


def read_trophic_groups(fpath, taxonomy):
    trophic_groups = pd.read_excel(fpath)
    taxa_to_trophic_group = dict()

    for row in trophic_groups.to_dict(orient="records"):
        if isinstance(row["Trophic group"], float) and isnan(row["Trophic group"]):
            continue
        else:
            for i in reversed(range(1, 6)):
                if not (
                    isinstance(row[f"Taxa_{i}"], float) and isnan(row[f"Taxa_{i}"])
                ):
                    taxa_to_trophic_group[row[f"Taxa_{i}"].strip()] = (
                        row["Trophic group"].strip().replace("parasitic", "parasite")
                    )
    group_to_taxa = defaultdict(set)
    for taxa, group in taxa_to_trophic_group.items():
        group_to_taxa[group].add(taxa)
    otu_to_group = dict()
    for otu_id, lineage in taxonomy.items():
        for trophic_group, taxons in group_to_taxa.items():
            if set(lineage) & taxons:
                otu_to_group[otu_id] = trophic_group
    return otu_to_group


def read_experimental_interactions(interactions_fpath: str) -> dict:
    interactions = dict()
    with open(interactions_fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for head, tail, _, habitat, interaction in reader:
            head_lineage = tuple(json.loads(head))
            tail_lineage = tuple(json.loads(tail))
            interactions[tuple(sorted([head_lineage, tail_lineage]))] = {
                "interaction": interaction,
                "habitat": habitat,
            }
    return interactions


def read_predicted_interactions(interactions_fpath: str) -> dict:
    interactions = dict()
    with open(interactions_fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for head, tail, _, _, sign in reader:
            interactions[tuple(sorted([head, tail]))] = {"sign": sign}
    return interactions


def read_taxonomy(tax_fpath: str) -> dict:
    taxonomy = dict()
    with open(tax_fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for otu_id, *lineage in reader:
            taxonomy[otu_id] = tuple(lineage)
    return taxonomy


def get_prop_known_interactions(graph, known_interactions, taxonomy):
    n = 0
    for head, tail in graph.edges():
        if head in taxonomy and tail in taxonomy:
            head_lineage = taxonomy[head]
            tail_lineage = taxonomy[tail]
            if tuple(sorted([head_lineage, tail_lineage])) in known_interactions:
                n += 1
    return n / len(known_interactions)


def get_prop_predicted_interactions(graph, predicted_interactions):
    n = 0
    for head, tail in graph.edges():
        if tuple(sorted([head, tail])) in predicted_interactions:
            n += 1
    return n / len(predicted_interactions)


def get_n_signs(graph: nx.Graph, sign: str) -> int:
    n = 0
    for _, _, attrs in graph.edges(data=True):
        if attrs["sign"] == sign:
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
@return_if_except(ValueError, np.inf)
def get_radius(graph: nx.Graph) -> Union[int, float]:
    return nx.radius(graph)


@return_if_except(nx.NetworkXError, np.inf)
@return_if_except(ValueError, np.inf)
def get_diameter(graph: nx.Graph) -> Union[int, float]:
    return nx.diameter(graph)


@return_if_except(nx.NetworkXError, np.nan)
def get_avg_path_length(graph: nx.Graph) -> Union[int, float]:
    return nx.average_shortest_path_length(graph)


@return_if_except(ZeroDivisionError, np.nan)
def get_avg_clustering(graph: nx.Graph) -> Union[int, float]:
    return nx.average_clustering(graph)


def get_trophic_group_relations(
    graph: nx.Graph, trophic_groups: Dict[str, str], trophic_pair: Tuple[str, str]
) -> int:
    n = 0
    for head, tail in graph.edges():
        if head in trophic_groups and tail in trophic_groups:
            groups = [trophic_groups[head], trophic_groups[tail]]
            if tuple(sorted(groups)) == trophic_pair:
                n += 1
    return n


METRIC_TO_FUN: Dict[str, Callable[[nx.Graph], Union[float, int]]] = {
    "Number of nodes": lambda graph: graph.number_of_nodes(),
    "Number of edges": lambda graph: graph.number_of_edges(),
    "Number of positive interactions": partial(get_n_signs, sign="+"),
    "Number of negative interactions": partial(get_n_signs, sign="-"),
    "Number of connected components": lambda graph: nx.number_connected_components(
        graph
    ),
    "Mean network centrality": lambda graph: np.mean(
        list(nx.degree_centrality(graph).values())
    ),
    "Mean network betweeness": lambda graph: np.mean(
        list(nx.betweenness_centrality(graph).values())
    ),
    "Network radius": get_radius,
    "Network diameter": get_diameter,
    "Average clustering coefficient": get_avg_clustering,
    "Average path length": get_avg_path_length,
}

HEADER = list(METRIC_TO_FUN.keys())


def extend_metrics(
    known_interactions, predicted_interactions, taxonomy, trophic_groups
):
    METRIC_TO_FUN["Discovered known interactions"] = partial(
        get_prop_known_interactions,
        known_interactions=known_interactions,
        taxonomy=taxonomy,
    )
    HEADER.append("Discovered known interactions")
    METRIC_TO_FUN["Discovered Lima-Mendez interactions"] = partial(
        get_prop_predicted_interactions, predicted_interactions=predicted_interactions
    )
    HEADER.append("Discovered Lima-Mendez interactions")
    all_groups = sorted(set(trophic_groups.values()))
    for group1, group2 in combinations(all_groups, 2):
        METRIC_TO_FUN[f"{group1}-{group2} interactions"] = partial(
            get_trophic_group_relations,
            trophic_groups=trophic_groups,
            trophic_pair=tuple(sorted([group1, group2])),
        )
        HEADER.append(f"{group1}-{group2} interactions")


METHOD_CAPITALIZATION = {
    "fastspar_results": "FastSpar (SparCC)",
    "flashweave_results": "FlashWeave",
    "phyloseq_results": "phyloseq",
    "conet_results": "Custom CoNet",
}


def prepare_stats(graph: nx.Graph):
    return {name: fun(graph) for name, fun in METRIC_TO_FUN.items()}


def read_graph(fpath: str) -> nx.Graph:
    graph = nx.Graph()
    with open(fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for node1, node2, *attrs in reader:
            if node1 != "leftover_vector" and node2 != "leftover_vector":
                graph.add_edge(node1, node2, sign=attrs[1], weight=attrs[0])
    return graph


def process(fpath: str):
    graph = read_graph(fpath)
    return fpath, prepare_stats(graph)


def write_interactions(
    graph: nx.Graph,
    graph_name: str,
    handle: OpenPyXLWriter,
    taxonomy: Dict[str, Tuple[str, ...]],
    known_interactions,
    predicted_interactions,
    trophic_groups,
):
    handle.writerow(
        [
            "otu 1",
            "otu 2",
            "weight",
            "sign",
            "in PIDA",
            "PIDA interaction",
            "PIDA habitat",
            "in Lima-Mendez",
            "Lima-Mendez sign",
            "trophic group 1",
            "trophic group 2",
        ]
        + list(chain(*[[f"{rank} {i}" for rank in RANKS] for i in range(1, 3)])),
        graph_name,
    )
    for head, tail, attrs in graph.edges(data=True):
        head, tail = sorted([head, tail])
        head_taxonomy = taxonomy[head]
        tail_taxonomy = taxonomy[tail]
        if tuple(sorted([head_taxonomy, tail_taxonomy])) in known_interactions:
            known = True
            habitat = known_interactions[tuple(sorted([head_taxonomy, tail_taxonomy]))][
                "habitat"
            ]
            interaction = known_interactions[
                tuple(sorted([head_taxonomy, tail_taxonomy]))
            ]["interaction"]
        else:
            known = False
            habitat = ""
            interaction = ""
        if (head, tail) in predicted_interactions:
            predicted = True
            predicted_sign = predicted_interactions[(head, tail)]["sign"]
        else:
            predicted = False
            predicted_sign = ""
        row = [
            head,
            tail,
            attrs["weight"],
            attrs["sign"],
            known,
            habitat,
            interaction,
            predicted,
            predicted_sign,
            trophic_groups.get(head, ""),
            trophic_groups.get(tail, ""),
            *(
                list(head_taxonomy)
                + ["" for _ in range(len(RANKS) - len(head_taxonomy))]
            ),
            *(
                list(tail_taxonomy)
                + ["" for _ in range(len(RANKS) - len(head_taxonomy))]
            ),
        ]
        handle.writerow(row, graph_name)


if __name__ == "__main__":
    proc_delayed = delayed(process)
    executor = Parallel(n_jobs=snakemake.threads, verbose=5)
    taxonomy = read_taxonomy(snakemake.input[0])
    known_interactions = read_experimental_interactions(snakemake.input[1])
    predicted_interactions = read_predicted_interactions(snakemake.input[2])
    trophic_groups = read_trophic_groups(snakemake.input[3], taxonomy)
    extend_metrics(known_interactions, predicted_interactions, taxonomy, trophic_groups)
    with OpenPyXLWriter(snakemake.output[0], "Network statistics") as handle:
        handle.writerow(["network"] + HEADER)
        tasks = (proc_delayed(fpath) for fpath in snakemake.input[4:])
        for fpath, stats in executor(tasks):
            if "consensus_network" in fpath:
                network_name = "Consensus network"
            else:
                head, tail = os.path.split(fpath)
                method_name = METHOD_CAPITALIZATION[os.path.split(head)[1]]
                network_name = method_name + " network"
            handle.writerow(
                [network_name] + [stats[statistic] for statistic in HEADER],
                "Network statistics",
            )
            write_interactions(
                read_graph(fpath),
                network_name,
                handle,
                taxonomy,
                known_interactions,
                predicted_interactions,
                trophic_groups,
            )