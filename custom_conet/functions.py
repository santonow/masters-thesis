import os
from operator import itemgetter
from typing import Optional, Tuple, Union, Iterable, Any, Set
from random import sample

import numpy as np
from numba import njit
from scipy.stats import chi2


MIN_DOUBLE = np.finfo(np.float64).min
SMALL_FLOAT = 0.0000001


def write_matrix(dirname: str, matrix: np.ndarray, name: str):
    np.save(os.path.join(dirname, name), matrix)


@njit
def std(arr, arr_mean) -> float:
    diff = arr - arr_mean
    return np.sqrt(np.sum(diff**2) / (arr.shape[0] - 1))


@njit
def pearson(x: np.ndarray, y: np.ndarray) -> float:
    n = x.shape[0]
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    std_x = std(x, x_mean)
    std_y = std(y, y_mean)
    return np.sum(((x - x_mean) / std_x) * ((y - y_mean) / std_y)) / (n - 1)


@njit
def get_ranks(arr: np.ndarray):
    """https://stackoverflow.com/a/5284703/12931685"""
    temp = arr.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(arr))
    return ranks + 1


@njit
def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(get_ranks(x), get_ranks(y))


@njit
def braycurtis(x: np.ndarray, y: np.ndarray) -> float:
    dist = 0.0
    total = 0.0
    for _x, _y in zip(x, y):
        dist += min(_x, _y)
        total += _x + _y
    return 1.0 - (2.0 * dist) / total


@njit
def kullback_leibler(x: np.ndarray, y: np.ndarray) -> float:
    # symmetric variant
    new_x = x + 1e-09
    new_y = y + 1e-09
    new_x /= sum(new_x)
    new_y /= sum(new_y)
    return float(np.sum(new_x * np.log(new_x / new_y))) + float(
        np.sum(new_y * np.log(new_y / new_x))
    )


@njit
def compute_correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if method == "pearson":
        return pearson(x, y)
    if method == "spearman":
        return spearman(x, y)
    if method == "bray-curtis":
        return braycurtis(x, y)
    if method == "kullback-leibler":
        return kullback_leibler(x, y)
    return np.nan


@njit
def compute_correlations(
    matrix: np.ndarray,
    method: str,
    xs: Optional[np.ndarray] = None,
    ys: Optional[np.ndarray] = None,
) -> np.ndarray:
    correlations = np.zeros((matrix.shape[1], matrix.shape[1]))
    if xs is not None and ys is not None:
        for i, j in zip(xs, ys):
            correlations[i, j] = compute_correlation(matrix[:, i], matrix[:, j], method)
            correlations[j, i] = correlations[i, j]
    else:
        for i in range(matrix.shape[1]):
            for j in range(matrix.shape[1]):
                if i < j and i != j:
                    correlations[i, j] = compute_correlation(
                        matrix[:, i], matrix[:, j], method
                    )
                    correlations[j, i] = correlations[i, j]
    return correlations


@njit
def bootstrap(matrix: np.ndarray) -> np.ndarray:
    bs_matrix = np.zeros(matrix.shape)
    for i in range(matrix.shape[0]):
        bs_matrix[i, :] = matrix[np.random.choice(np.arange(0, matrix.shape[0])), :]
    for i in range(bs_matrix.shape[1]):
        if np.all(bs_matrix[:, i] == bs_matrix[0, i]):
            return bootstrap(matrix)
    return bs_matrix


@njit
def update_mean_variance(
    new_value: Union[np.ndarray, float],
    means: np.ndarray,
    variances: np.ndarray,
    count: int,
    indices: Optional[Tuple[int, int]] = None,
) -> None:
    """Welford's method for computing variance.

    This is an update step. First means and variances matrices have to be initialized to zeros.
    Inspired by: https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/.
    """
    if indices is None:
        prev_means = means.copy()
        means += (new_value - means) / count
        variances += (new_value - means) * (new_value - prev_means)
    else:
        prev_means = means[indices]
        means[indices] += (new_value - means[indices]) / count
        variances[indices] += (new_value - means[indices]) * (new_value - prev_means)


@njit
def compute_p_value(
    permutation_means: np.ndarray,
    permutation_variances: np.ndarray,
    bootstrap_means: np.ndarray,
    bootstrap_variances: np.ndarray,
) -> np.ndarray:
    z_scores = (bootstrap_means - permutation_means) / np.sqrt(
        0.5 * (permutation_variances + bootstrap_variances)
    )
    p_values = 2 * (np.exp(-((-np.abs(z_scores)) ** 2)) / 2) / np.sqrt(2 * np.pi)
    return p_values


@njit
def normalize(matrix: np.ndarray):
    for i in range(matrix.shape[0]):
        s = np.sum(matrix[i])
        for j in range(matrix.shape[1]):
            matrix[i, j] = matrix[i, j] / s


@njit
def permute(
    matrix: np.ndarray,
    permutation_means: np.ndarray,
    permutation_variances: np.ndarray,
    method: str,
    iteration: int,
    xs: Optional[np.ndarray] = None,
    ys: Optional[np.ndarray] = None,
    renorm: bool = True,
):
    diffs_matrix = np.zeros(matrix.shape)
    if renorm:
        for i in range(matrix.shape[0]):
            s = matrix[i].sum()
            for j in range(matrix.shape[1]):
                diffs_matrix[i, j] = s - matrix[i, j]

    if xs is not None and ys is not None:
        for i, j in zip(xs, ys):
            col1 = np.random.permutation(matrix[:, i])
            col2 = np.random.permutation(matrix[:, j])
            if renorm:
                col1 /= diffs_matrix[:, i] + col1
                col2 /= diffs_matrix[:, j] + col2
            corr = compute_correlation(col1, col2, method)
            update_mean_variance(
                corr, permutation_means, permutation_variances, iteration, (i, j)
            )
            permutation_means[j, i] = permutation_means[i, j]
            permutation_variances[j, i] = permutation_variances[i, j]
    else:
        for i in range(matrix.shape[1]):
            for j in range(matrix.shape[1]):
                if i < j and i != j:
                    col1 = np.random.permutation(matrix[:, i])
                    col2 = np.random.permutation(matrix[:, j])
                    if renorm:
                        col1 /= diffs_matrix[:, i] + col1
                        col2 /= diffs_matrix[:, j] + col2
                    corr = compute_correlation(col1, col2, method)
                    update_mean_variance(
                        corr,
                        permutation_means,
                        permutation_variances,
                        iteration,
                        (i, j),
                    )
                    permutation_means[j, i] = permutation_means[i, j]
                    permutation_variances[j, i] = permutation_variances[i, j]


def reboot(
    matrix: np.ndarray,
    n_iter: int,
    method: str,
    indices: Set[Tuple[int, int]],
    renorm: bool = False,
    samples_for_ci: int = 100,
):
    """Perform ReBoot procedure."""
    # for determining a confidence interval
    samples = []
    xs = np.array([x[0] for x in sorted(indices)], dtype=np.int64)
    ys = np.array([x[1] for x in sorted(indices)], dtype=np.int64)

    # bootstrap
    bs_means = np.zeros((matrix.shape[1], matrix.shape[1]))
    bs_variances = np.zeros((matrix.shape[1], matrix.shape[1]))
    for i in range(n_iter):
        bs_matrix = bootstrap(matrix)
        bs_correlations = compute_correlations(bs_matrix, method, xs, ys)
        samples.extend(
            sample(
                [x[0] for x in matrix_iter(bs_correlations, indices)], samples_for_ci
            )
        )
        update_mean_variance(bs_correlations, bs_means, bs_variances, i + 1)

    # permutations
    permutation_means = np.zeros((matrix.shape[1], matrix.shape[1]))
    permutation_variances = np.zeros((matrix.shape[1], matrix.shape[1]))
    if method == "kullback-leibler":
        renorm = False
    for i in range(n_iter):
        permute(
            matrix,
            permutation_means,
            permutation_variances,
            method,
            i + 1,
            xs,
            ys,
            renorm,
        )
    return (
        permutation_means,
        permutation_variances,
        bs_means,
        bs_variances,
        n_iter,
        samples,
    )


def compute_pvals(matrix: np.ndarray, n_iter: int, method: str, indices):
    (
        permutation_means,
        permutation_variances,
        bs_means,
        bs_variances,
        count,
        bs_samples,
    ) = reboot(matrix, n_iter, method, indices=indices)
    permutation_variances /= count
    bs_variances = bs_variances / count
    pvals = compute_p_value(
        permutation_means, permutation_variances, bs_means, bs_variances
    )

    return pvals, bs_samples


@njit
def update_mean_variance_parallel(
    prev_n: int,
    prev_means: np.ndarray,
    prev_vars: np.ndarray,
    new_n: int,
    new_means: np.ndarray,
    new_vars: np.ndarray,
):
    n = prev_n + new_n
    delta = new_means - prev_means
    _vars = prev_vars + new_vars + (delta**2) * prev_n * new_n / n
    means = (prev_n * prev_means + new_n * new_means) / n
    return means, _vars


@njit
def clr(matrix: np.ndarray):
    matrix = matrix + SMALL_FLOAT
    geom_mean = np.zeros(matrix.shape[0])
    for i in range(matrix.shape[0]):
        geom_mean[i] = np.exp(np.mean(np.log(matrix[i, :])))
    for i in range(matrix.shape[0]):
        matrix[i, :] /= geom_mean[i]
    return np.log(matrix)


def matrix_iter(
    matrix: np.ndarray, indices: Optional[Set[Tuple[int, int]]] = None
) -> Iterable[Tuple[Any, int, int]]:
    if indices is None:
        all_indices = True
    else:
        all_indices = False
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if (
                i != j
                and i < j
                and (all_indices or (i, j) in indices or (j, i) in indices)
            ):
                yield matrix[i, j], i, j


@njit
def fisher_merge(
    pvals_1: np.ndarray, pvals_2: np.ndarray, xs: np.ndarray, ys: np.ndarray
) -> np.ndarray:
    merged_pvals = np.empty_like(pvals_1)
    if xs is not None and ys is not None:
        for i, j in zip(xs, ys):
            merged_pval = 0
            for m in [pvals_1, pvals_2]:
                pval = m[i, j]
                if pval == 0.0:
                    pval = MIN_DOUBLE
                merged_pval += np.log(pval)
            merged_pvals[i, j] = merged_pval
            merged_pvals[j, i] = merged_pval
    else:
        for i in range(pvals_1.shape[0]):
            for j in range(pvals_1.shape[0]):
                if i < j and i != j:
                    merged_pval = 0
                    for m in [pvals_1, pvals_2]:
                        pval = m[i, j]
                        if pval == 0.0:
                            pval = MIN_DOUBLE
                        merged_pval += np.log(pval)
                    merged_pvals[i, j] = merged_pval
                    merged_pvals[j, i] = merged_pval
    return (-2.0) * merged_pvals


def chi_square_p_vals(m: np.ndarray, dof) -> np.ndarray:
    result = np.empty_like(m)
    for i in range(m.shape[0]):
        for j in range(m.shape[0]):
            if i < j and i != j:
                result[i, j] = chi2.sf(m[i, j], dof)
                result[j, i] = result[i, j]
    return result


def merge_p_values(pvals_1, pvals_2, indices):
    """Brown's method of merging correlations.

    Adapted from CONet code.
    """
    # first, get an approx correction factor (which is a modified covariance of two p-value sets) and dof
    measure_number = 2
    pvals_1_vec = np.array(
        [pval for pval, _, _ in matrix_iter(pvals_1, indices)], dtype=np.float64
    )
    pvals_2_vec = np.array(
        [pval for pval, _, _ in matrix_iter(pvals_2, indices)], dtype=np.float64
    )
    corrcoeff = compute_correlation(pvals_1_vec, pvals_2_vec, method="pearson")

    # Brown's approx formula
    if corrcoeff > 0:
        variance = corrcoeff * (3.25 + 0.75 * corrcoeff)
    else:
        variance = corrcoeff * (3.27 + 0.71 * corrcoeff)

    variance = 2 * variance
    variance = variance + measure_number / 4

    dof = (2 * ((measure_number * 2) ** 2)) / variance
    correction_factor = variance / (2.0 * 2 * measure_number)
    # compute corrected p-values

    xs = np.array([x[0] for x in sorted(indices)], dtype=np.int64)
    ys = np.array([x[1] for x in sorted(indices)], dtype=np.int64)

    chi_square = fisher_merge(pvals_1, pvals_2, xs, ys) / correction_factor
    return chi_square_p_vals(chi_square, dof)


def benjamini_hochberg(matrix: np.ndarray, indices):
    corrected = np.empty_like(matrix)
    n_pvals = matrix.shape[0] * matrix.shape[1]
    for i, (pval, row, col) in enumerate(
        sorted(matrix_iter(matrix, indices), key=itemgetter(0))
    ):
        corrected[row, col] = pval * n_pvals / (i + 1)
        corrected[col, row] = corrected[row, col]
    return corrected
