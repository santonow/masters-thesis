import json
from multiprocessing.spawn import prepare
import os
from collections import defaultdict, Counter
from functools import partial
from itertools import combinations, chain
import sys
from typing import Dict, Callable, Union, Any, List, Optional, Tuple
import csv
from math import isnan

import networkx as nx
import numpy as np
import openpyxl
from openpyxl.utils import get_column_letter
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
        i = 1
        while True:
            yield f"{get_column_letter(i)}{self._current_pos}"
            i += 1

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


def prepare_pairwise_stats(
    method_graphs: dict[str, nx.Graph],
    handle: OpenPyXLWriter,
    taxonomy: dict[str, tuple[str, ...]],
):
    def common_taxon_edges(c1: Counter, c2: Counter) -> int:
        common = 0
        for k in set(c1) & set(c2):
            common += min([c1[k], c2[k]])
        return common

    edges_OTU = {
        method: set(tuple(sorted([head, tail])) for head, tail in graph.edges())
        for method, graph in method_graphs.items()
    }
    edges_taxon = {
        method: Counter(
            tuple(sorted([taxonomy[head], taxonomy[tail]]))
            for head, tail in graph.edges()
        )
        for method, graph in method_graphs.items()
    }

    common_OTUs = {
        tuple(sorted([method1, method2])): len(edges_OTU[method1] & edges_OTU[method2])
        for method1, method2 in combinations(method_graphs, 2)
    }
    common_taxons = {
        tuple(sorted([method1, method2])): common_taxon_edges(
            edges_taxon[method1], edges_taxon[method2]
        )
        for method1, method2 in combinations(method_graphs, 2)
    }

    methods = sorted(method_graphs)
    handle.writerow([""] + methods, sheet_name="Method comparisons")
    for level, common_edges in [
        ("OTU level", common_OTUs),
        ("taxon level", common_taxons),
    ]:
        for method in methods:
            row = [f"{method} ({level})"]
            for other_method in methods:
                if method == other_method:
                    row.append("")
                else:
                    row.append(common_edges[tuple(sorted([method, other_method]))])

            handle.writerow(row, sheet_name="Method comparisons")


def read_trophic_groups(fpath, taxonomy):
    trophic_groups = pd.read_excel(fpath)
    taxa_to_trophic_group = dict()

    for row in trophic_groups.to_dict(orient="records"):
        if isinstance(row["Trophic group"], float):
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


def read_predicted_interactions(interactions_fpath: str) -> dict[str, dict[str, str]]:
    interactions = dict()
    with open(interactions_fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for head, tail, _, _, sign in reader:
            interactions[tuple(sorted([head, tail]))] = {"sign": sign}
    return interactions


def read_taxonomy(tax_fpath: str) -> dict[str, tuple[str, ...]]:
    taxonomy = dict()
    with open(tax_fpath) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for otu_id, *lineage in reader:
            taxonomy[otu_id] = tuple(lineage)
    return taxonomy


def get_prop_known_interactions(
    graph: nx.Graph,
    known_interactions,
    taxonomy: dict[str, tuple[str, ...]],
    count: bool = False,
    prop_of_network_edges: bool = False,
    kind="all",
) -> float:
    n = 0
    results = Counter()
    all_tax_relations = set()
    for head, tail in graph.edges():
        if head in taxonomy and tail in taxonomy:
            head_lineage = taxonomy[head]
            tail_lineage = taxonomy[tail]
            tax_relation = tuple(sorted([head_lineage, tail_lineage]))
            all_tax_relations.add(tax_relation)
            if tax_relation in known_interactions:
                data = known_interactions[tax_relation]
                results["all"] += 1
                results[data["interaction"]] += 1
                n += 1
    if count:
        return results[kind]
    else:
        if prop_of_network_edges:
            return results[kind] / len(all_tax_relations)
        else:
            return results[kind] / len(known_interactions)


def get_prop_predicted_interactions(
    graph: nx.Graph,
    predicted_interactions: dict[str, dict[str, str]],
    taxonomy: Optional[dict[str, tuple[str, ...]]],
    count: bool = False,
    prop_of_network_edges: bool = False,
) -> float:
    if taxonomy is not None:
        inferred_interactions = {
            tuple(sorted([taxonomy[head], taxonomy[tail]]))
            for head, tail in graph.edges()
        }
        interactions = {
            tuple(sorted([taxonomy[head], taxonomy[tail]]))
            for head, tail in predicted_interactions
        }
    else:
        inferred_interactions = {
            tuple(sorted([head, tail])) for head, tail in graph.edges()
        }
        interactions = set(predicted_interactions)
    if count:
        return len(inferred_interactions & interactions)
    else:
        if prop_of_network_edges:
            return len(inferred_interactions & interactions) / len(
                inferred_interactions
            )
        else:
            return len(inferred_interactions & interactions) / len(interactions)


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
    "Mean network degree centrality": lambda graph: np.mean(
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
        count=True,
    )
    HEADER.append("Discovered known interactions")

    METRIC_TO_FUN["Proportion of known interactions discovered"] = partial(
        get_prop_known_interactions,
        known_interactions=known_interactions,
        taxonomy=taxonomy,
    )
    HEADER.append("Proportion of known interactions discovered")

    METRIC_TO_FUN["Proportion of edges that appear as known interactions"] = partial(
        get_prop_known_interactions,
        known_interactions=known_interactions,
        taxonomy=taxonomy,
        prop_of_network_edges=True,
    )
    HEADER.append("Proportion of edges that appear as known interactions")

    all_interactions = {elem["interaction"] for elem in known_interactions.values()}
    for kind in sorted(all_interactions):
        METRIC_TO_FUN[f"Discovered known interactions ({kind})"] = partial(
            get_prop_known_interactions,
            known_interactions=known_interactions,
            taxonomy=taxonomy,
            count=True,
            kind=kind,
        )
        HEADER.append(f"Discovered known interactions ({kind})")

    METRIC_TO_FUN["Discovered Lima-Mendez interactions (OTU level)"] = partial(
        get_prop_predicted_interactions,
        predicted_interactions=predicted_interactions,
        taxonomy=None,
        count=True,
    )
    HEADER.append("Discovered Lima-Mendez interactions (OTU level)")

    METRIC_TO_FUN[
        "Proportion of Lima-Mendez interactions (OTU level) discovered"
    ] = partial(
        get_prop_predicted_interactions,
        predicted_interactions=predicted_interactions,
        taxonomy=None,
    )
    HEADER.append("Proportion of Lima-Mendez interactions (OTU level) discovered")

    METRIC_TO_FUN[
        "Proportion of edges that appear as Lima-Mendez interactions (OTU level)"
    ] = partial(
        get_prop_predicted_interactions,
        predicted_interactions=predicted_interactions,
        taxonomy=None,
        prop_of_network_edges=True,
    )
    HEADER.append(
        "Proportion of edges that appear as Lima-Mendez interactions (OTU level)"
    )

    METRIC_TO_FUN["Discovered Lima-Mendez interactions (taxonomy level)"] = partial(
        get_prop_predicted_interactions,
        predicted_interactions=predicted_interactions,
        taxonomy=taxonomy,
        count=True,
    )
    HEADER.append("Discovered Lima-Mendez interactions (taxonomy level)")

    METRIC_TO_FUN[
        "Proportion of Lima-Mendez interactions (taxonomy level) discovered"
    ] = partial(
        get_prop_predicted_interactions,
        predicted_interactions=predicted_interactions,
        taxonomy=taxonomy,
    )
    HEADER.append("Proportion of Lima-Mendez interactions (taxonomy level) discovered")

    METRIC_TO_FUN[
        "Proportion of edges that appear as Lima-Mendez interactions (taxonomy level)"
    ] = partial(
        get_prop_predicted_interactions,
        predicted_interactions=predicted_interactions,
        taxonomy=taxonomy,
        prop_of_network_edges=True,
    )
    HEADER.append(
        "Proportion of edges that appear as Lima-Mendez interactions (taxonomy level)"
    )

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
    return fpath, prepare_stats(graph), graph


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
        key = tuple(sorted([head_taxonomy, tail_taxonomy]))
        if key in known_interactions:
            known = True
            habitat = known_interactions[key]["habitat"]
            interaction = known_interactions[key]["interaction"]
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
                + ["" for _ in range(len(RANKS) - len(tail_taxonomy))]
            ),
        ]
        handle.writerow(row, graph_name)


if __name__ == "__main__":
    (
        tax_path,
        known_inter_path,
        pred_path,
        trophic_groups_path,
        n_threads,
        output_fpath,
        *network_fpaths,
    ) = sys.argv[1:]
    proc_delayed = delayed(process)
    executor = Parallel(n_jobs=int(n_threads), verbose=5)
    taxonomy = read_taxonomy(tax_path)
    known_interactions = read_experimental_interactions(known_inter_path)
    predicted_interactions = read_predicted_interactions(pred_path)
    trophic_groups = read_trophic_groups(trophic_groups_path, taxonomy)
    extend_metrics(known_interactions, predicted_interactions, taxonomy, trophic_groups)
    tasks = (proc_delayed(fpath) for fpath in network_fpaths)
    graphs_by_method = {}
    stats_by_method = {}
    for fpath, stats, graph in executor(tasks):
        if "consensus_network" in fpath:
            network_name = "Consensus network"
        else:
            head, tail = os.path.split(fpath)
            method_name = METHOD_CAPITALIZATION[os.path.split(head)[1]]
            network_name = method_name + " network"
        graphs_by_method[network_name] = graph
        stats_by_method[network_name] = stats

    with OpenPyXLWriter(output_fpath, "Network statistics") as handle:
        handle.writerow(["network"] + HEADER)
        for network_name, stats in stats_by_method.items():
            handle.writerow(
                [network_name] + [stats[statistic] for statistic in HEADER],
                "Network statistics",
            )
        prepare_pairwise_stats(graphs_by_method, handle, taxonomy)
        for network_name, graph in graphs_by_method.items():
            write_interactions(
                graph,
                network_name,
                handle,
                taxonomy,
                known_interactions,
                predicted_interactions,
                trophic_groups,
            )
