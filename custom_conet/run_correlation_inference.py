import csv
import sys
from collections import Counter
from typing import Tuple, Set

from joblib import Parallel
from joblib import delayed as delayed
import numpy as np
import networkx as nx

from occurence_table import OTUTable
from functions import (
    compute_correlations,
    compute_p_value,
    reboot,
    compute_pvals,
    update_mean_variance_parallel,
    normalize,
    benjamini_hochberg,
    merge_p_values,
    clr,
    matrix_iter,
)


class timer:
    def __init__(self, message):
        self.message = message

    def __enter__(self):
        print(self.message)
        self.start_time = time.time()

    def __exit__(self, exc_type, exc_val, exc_tb):
        print(f"Done in {time.time() - self.start_time:.4f}s.")
        print()


class CoNet:
    def __init__(
        self,
        otu_table: OTUTable,
        n_iter: int = 1000,
        methods: Tuple[str] = ("spearman", "kullback-leibler"),
        p_value_threshold: float = 0.05,
        n_jobs: int = 4,
        bh: bool = True,
        confidence_interval: float = 0.95,
        n_initial_edges: int = 10000000,
    ):
        self.otu_table = OTUTable
        self.indice_to_id = {
            i: _id for i, _id in enumerate(otu_table.ids(axis="observation"))
        }
        self.matrix = otu_table.matrix.T
        self.n_iter = n_iter
        self.methods = methods
        self.p_value_threshold = p_value_threshold
        self.n_jobs = n_jobs
        self.bh = bh
        self.confidence_interval = confidence_interval
        self.n_initial_edges = n_initial_edges
        self.matrices_by_method = self.prepare_matrices()
        by_method, common_edges = self.initial_correlation()
        self.correlations_by_method = by_method
        self.common_edges = common_edges
        self.results = None

    def prepare_matrices(self):
        matrices = dict()
        for method in self.methods:
            if method != "kullback-leibler":
                with timer(f"Computing clr transformation for method {method}..."):
                    matrix = clr(self.matrix)
            else:
                matrix = self.matrix.copy()
            if method == "kullback-leibler":
                with timer("Normalizing..."):
                    normalize(matrix)
            matrices[method] = matrix
        return matrices

    def initial_correlation(self):
        all_edges = []
        correlations_by_method = dict()
        for method, matrix in self.matrices_by_method.items():
            with timer(f"Computing correlations with method {method}..."):
                correlations = compute_correlations(matrix, method)
            correlations_by_method[method] = correlations
            edges, left, right = self.select_edges(
                correlations, method, self.n_initial_edges
            )
            all_edges.append(edges)
            print(
                f"Selecting {self.n_initial_edges} initial edges, metric value < {left} and > {right}"
            )
        common_edges = set.intersection(*all_edges)
        print(f"Selecting {len(common_edges)} common initial edges.")
        return correlations_by_method, common_edges

    def run(self):
        results = dict()
        for method, matrix in self.matrices_by_method.items():
            with timer(
                f"Computing p-values for method {method} using {self.n_jobs} processes..."
            ):
                if self.n_jobs == 1:
                    pvalues, bs_samples = compute_pvals(
                        matrix, self.n_iter, method, self.common_edges
                    )
                else:
                    pvalues, bs_samples = self.compute_pvals_parallel(
                        matrix=matrix, method=method
                    )
                left_interval, right_interval = np.percentile(
                    np.array(bs_samples),
                    [
                        100 * (1 - self.confidence_interval) / 2,
                        100 * (1 - (1 - self.confidence_interval) / 2),
                    ],
                )
                pvalues[pvalues > 0.5] = (
                    1 - pvalues[pvalues > 0.5]
                )  # we are doing one-sided test
                results[method] = {
                    "correlations": self.correlations_by_method[method],
                    "pvalues": pvalues,
                    "confidence_interval": (left_interval, right_interval),
                }
        if len(results) > 1:
            with timer("Merging p values..."):
                pvals = merge_p_values(
                    *[result["pvalues"] for result in results.values()],
                    self.common_edges,
                )
        else:
            pvals = next(iter(results.values()))["pvalues"]
        if self.bh:
            with timer("Performing Benjamini-Hochberg correction..."):
                pvals = benjamini_hochberg(pvals, self.common_edges)

        self.results = {"by_method": results, "pvals": pvals}
        rejected_edges = Counter()
        network = nx.Graph()
        for method, result in results.items():
            correlations = result["correlations"]
            left_interval, right_interval = result["confidence_interval"]
            method_pvalues = result["pvalues"]
            for i in range(correlations.shape[0]):
                for j in range(correlations.shape[1]):
                    if (
                        i != j
                        and i < j
                        and ((i, j) in self.common_edges or (j, i) in self.common_edges)
                    ):
                        if pvals[i, j] < self.p_value_threshold and (
                            left_interval < correlations[i, j] < right_interval
                        ):
                            left_otu_id = self.indice_to_id[i]
                            right_otu_id = self.indice_to_id[j]
                            if not network.has_edge(left_otu_id, right_otu_id):
                                network.add_edge(
                                    left_otu_id,
                                    right_otu_id,
                                    pvalue=pvals[i, j],
                                )
                            network[left_otu_id][right_otu_id][
                                f"{method}_weight"
                            ] = correlations[i, j]
                            network[left_otu_id][right_otu_id][
                                f"{method}_pvalue"
                            ] = method_pvalues[i, j]
                        else:
                            rejected_edges[method] += 1
        print(
            f"Nodes in graph: {network.number_of_nodes()}, edges in network: {network.number_of_edges()}."
        )
        print("Rejected edges:", rejected_edges)

        return network

    @staticmethod
    def select_edges(
        matrix: np.ndarray, method: str, max_edges=100000
    ) -> Tuple[Set[Tuple[int, int]], int, int]:
        if method == "kullback-leibler":
            pos_edges = [
                (val, i, j)
                for val, i, j in sorted(matrix_iter(matrix), key=lambda x: x[0])[
                    : max_edges // 2
                ]
            ]
            neg_edges = [
                (val, i, j)
                for val, i, j in sorted(
                    matrix_iter(matrix), key=lambda x: x[0], reverse=True
                )[: max_edges // 2]
            ]
            left = max(x[0] for x in pos_edges)
            right = min(x[0] for x in neg_edges)
        else:
            pos_edges = [
                (val, i, j)
                for val, i, j in sorted(
                    matrix_iter(matrix), key=lambda x: x[0], reverse=True
                )[: max_edges // 2]
            ]
            neg_edges = [
                (val, i, j)
                for val, i, j in sorted(matrix_iter(matrix), key=lambda x: x[0])[
                    : max_edges // 2
                ]
            ]
            left = max(x[0] for x in neg_edges)
            right = min(x[0] for x in pos_edges)
        return (
            set([(i, j) for _, i, j in pos_edges] + [(i, j) for _, i, j in neg_edges]),
            left,
            right,
        )

    def compute_pvals_parallel(self, matrix: np.ndarray, method: str):
        def prepare_proc_iters(n_iter, n_jobs):
            iters = [n_iter // n_jobs for _ in range(n_jobs)]
            iters[-1] += n_iter - sum(iters)
            return iters

        executor = Parallel(n_jobs=self.n_jobs, verbose=10)
        delated_fun = delayed(reboot)
        tasks = (
            delated_fun(matrix, proc_iters, method, self.common_edges)
            for proc_iters in prepare_proc_iters(self.n_iter, self.n_jobs)
        )

        current_count = 0
        perm_means = np.zeros((matrix.shape[1], matrix.shape[1]))
        perm_variances = np.zeros((matrix.shape[1], matrix.shape[1]))
        bs_means = np.zeros((matrix.shape[1], matrix.shape[1]))
        bs_variances = np.zeros((matrix.shape[1], matrix.shape[1]))
        samples = []
        for result in executor(tasks):
            (
                curr_perm_means,
                curr_perm_variances,
                curr_bs_means,
                curr_bs_variances,
                count,
                proc_samples,
            ) = result
            perm_means, perm_variances = update_mean_variance_parallel(
                current_count,
                perm_means,
                perm_variances,
                count,
                curr_perm_means,
                curr_perm_variances,
            )
            bs_means, bs_variances = update_mean_variance_parallel(
                current_count,
                bs_means,
                bs_variances,
                count,
                curr_bs_means,
                curr_bs_variances,
            )
            current_count += count
            samples.extend(proc_samples)
        perm_variances /= current_count
        bs_variances /= current_count
        assert current_count == self.n_iter

        return (
            compute_p_value(perm_means, perm_variances, bs_means, bs_variances),
            samples,
        )


if __name__ == "__main__":
    import time

    inp, n_threads, output = sys.argv[1:4]
    otu_table = OTUTable.from_tsv(inp, sample_rows=False)
    start_time = time.time()
    methods = ["spearman", "kullback-leibler"]
    conet = CoNet(otu_table, n_jobs=int(n_threads), n_iter=1000)
    network = conet.run()
    print("Took", time.time() - start_time)
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["head", "tail", "merged pval"]
            + [f"{method}_weight" for method in methods]
            + [f"{method}_pvalue" for method in methods]
        )
        for head, tail, attrs in network.edges(data=True):
            if all(f"{method}_weight" in attrs for method in methods):
                writer.writerow(
                    [head, tail, attrs["pvalue"]]
                    + [attrs[f"{method}_weight"] for method in methods]
                    + [attrs[f"{method}_pvalue"] for method in methods]
                )
