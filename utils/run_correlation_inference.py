import csv
import sys
from typing import List, Union

from joblib import Parallel
from joblib import delayed as _delayed
import numpy as np
import networkx as nx

from occurence_table import OTUTable
from numba_functions import compute_correlations, compute_p_value, reboot, compute_pvals, \
    update_mean_variance_parallel, normalize, benjamini_hochberg, merge_p_values


def delayed(fun):
    return _delayed(fun)


@delayed
def permuter(matrix: np.ndarray, proc_iters: int, method: str):
    try:
        result = reboot(matrix, proc_iters, method)
    except Exception as e:
        print(matrix.shape, proc_iters, method)
        print(e)
        raise e
    return result


def compute_pvals_parallel(matrix: np.ndarray, n_iter: int, method: str, n_jobs: int):
    def prepare_proc_iters(n_iter, n_jobs):
        iters = [n_iter // n_jobs for _ in range(n_jobs)]
        iters[-1] += n_iter - sum(iters)
        return iters

    executor = Parallel(n_jobs=n_jobs, verbose=10)
    tasks = (permuter(matrix, proc_iters, method) for proc_iters in prepare_proc_iters(n_iter, n_jobs))

    current_count = 0
    perm_means = np.zeros((matrix.shape[1], matrix.shape[1]))
    perm_variances = np.zeros((matrix.shape[1], matrix.shape[1]))
    bs_means = np.zeros((matrix.shape[1], matrix.shape[1]))
    bs_variances = np.zeros((matrix.shape[1], matrix.shape[1]))
    samples = []
    for result in executor(tasks):
        curr_perm_means, curr_perm_variances, curr_bs_means, curr_bs_variances, count, proc_samples = result
        perm_means, perm_variances = update_mean_variance_parallel(
            current_count, perm_means, perm_variances,
            count, curr_perm_means, curr_perm_variances
        )
        bs_means, bs_variances = update_mean_variance_parallel(
            current_count, bs_means, bs_variances,
            count, curr_bs_means, curr_bs_variances
        )
        current_count += count
        samples.extend(proc_samples)
    perm_variances /= current_count
    bs_variances /= current_count
    assert current_count == n_iter

    return compute_p_value(perm_means, perm_variances, bs_means, bs_variances), samples


def infer_network(
    otu_table: OTUTable, n_iter: int, methods: Union[str, List[str]],
    p_value_threshold: float = 0.05, n_jobs: int = 4, bh: bool = True, confidence_interval=0.95
) -> nx.Graph:
    network = nx.Graph()
    indice_to_id = {i: _id for i, _id in enumerate(otu_table.ids(axis='observation'))}
    if isinstance(methods, str):
        methods = [methods]
    results = dict()
    for method in methods:
        matrix = otu_table.matrix.T
        if method != 'kullback-leibler':
            normalize(matrix)
        print(f'Computing correlations with method {method}...')
        correlations = compute_correlations(matrix, method)
        if n_jobs == 1:
            pvalues, bs_samples = compute_pvals(matrix, n_iter, method)
        else:
            pvalues, bs_samples = compute_pvals_parallel(matrix, n_iter, method, n_jobs)

        left_interval, right_interval = np.percentile(
            np.array(bs_samples),
            [100 * (1 - confidence_interval) / 2, 100 * (1 - (1 - confidence_interval) / 2)]
        )
        pvalues[pvalues > 0.5] = 1 - pvalues[pvalues > 0.5]  # we are doing one-sided test
        results[method] = {
            'correlations': correlations,
            'pvalues': pvalues,
            'confidence_interval': (left_interval, right_interval)
        }
    if len(results) > 1:
        pvals = merge_p_values(*[result['pvalues'] for result in results.values()])
    else:
        pvals = next(iter(results.values()))['pvalues']
    if bh:
        pvals = benjamini_hochberg(pvals)
    for method, result in results.items():
        correlations = result['correlations']
        left_interval, right_interval = result['confidence_interval']
        method_pvalues = result['pvalues']
        for i in range(correlations.shape[0]):
            for j in range(correlations.shape[1]):
                if i != j and i < j:
                    if (
                        pvals[i, j] < p_value_threshold
                        and (left_interval < correlations[i, j] < right_interval)
                    ):
                        left_otu_id = indice_to_id[i]
                        right_otu_id = indice_to_id[j]
                        if not network.has_edge(left_otu_id, right_otu_id):
                            network.add_edge(
                                left_otu_id, right_otu_id,
                                pvalue=pvals[i, j],
                            )
                        network[left_otu_id][right_otu_id][f'{method}_weight'] = correlations[i, j]
                        network[left_otu_id][right_otu_id][f'{method}_pvalue'] = method_pvalues[i, j]
    return network


if __name__ == '__main__':
    import time
    inp, n_threads, output = sys.argv[1:4]
    otu_table = OTUTable.from_tsv(inp, sample_rows=False)
    start_time = time.time()
    methods = ['kullback-leibler', 'spearman']
    network = infer_network(otu_table, 100, methods, 0.05, int(n_threads))
    print('Took', time.time() - start_time)
    with open(output, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(
            ['head', 'tail', 'merged pval']
            + [f'{method}_weight' for method in methods]
            + [f'{method}_pvalue' for method in methods]
        )
        for head, tail, attrs in network.edges(data=True):
            if all(f'{method}_weight' in attrs for method in methods):
                writer.writerow(
                    [head, tail, attrs['pvalue']]
                    + [attrs[f'{method}_weight'] for method in methods]
                    + [attrs[f'{method}_pvalue'] for method in methods]
                )
