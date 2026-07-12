"""Core replication-unit audit.

The functions operate on a long-form source-recipient table.  Each row is one
recipient measured near one perturbation source.  The implementation keeps the
point estimate fixed while comparing pair-naive and source-clustered uncertainty,
which is the central ablation reported in the manuscript.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class AuditResult:
    outcome: str
    estimate: float
    se_naive: float
    p_naive: float
    se_source: float
    p_source: float
    se_inflation: float
    n_pairs: int
    n_sources: int

    def to_dict(self) -> dict:
        return asdict(self)


def _design(frame: pd.DataFrame, treatment: str, covariates: Iterable[str]):
    columns = [treatment, *covariates]
    x = pd.get_dummies(frame.loc[:, columns], drop_first=True, dtype=float)
    if treatment not in x:
        raise ValueError(f"treatment column {treatment!r} was lost during encoding")
    x.insert(0, "const", 1.0)
    return x


def _ols(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    return beta, residual, xtx_inv


def _cluster_covariance(
    x: np.ndarray, residual: np.ndarray, xtx_inv: np.ndarray, groups: pd.Series
) -> tuple[np.ndarray, int]:
    codes, levels = pd.factorize(groups, sort=True)
    meat = np.zeros((x.shape[1], x.shape[1]))
    for code in range(len(levels)):
        mask = codes == code
        score = x[mask].T @ residual[mask]
        meat += np.outer(score, score)
    n, p, g = len(x), x.shape[1], len(levels)
    correction = (g / (g - 1)) * ((n - 1) / (n - p)) if g > 1 else 1.0
    return correction * xtx_inv @ meat @ xtx_inv, g


def compare_replication_units(
    frame: pd.DataFrame,
    outcome: str,
    treatment: str = "treat",
    source: str = "source_id",
    covariates: Iterable[str] = (),
) -> AuditResult:
    """Fit one OLS contrast with naive and source-clustered covariance."""
    needed = [outcome, treatment, source, *covariates]
    data = frame.dropna(subset=needed).copy()
    xdf = _design(data, treatment, covariates)
    x = xdf.to_numpy(float)
    y = data[outcome].astype(float).to_numpy()
    beta, residual, xtx_inv = _ols(y, x)
    n, p = x.shape
    sigma2 = float(residual @ residual / (n - p))
    cov_naive = sigma2 * xtx_inv
    cov_source, n_groups = _cluster_covariance(x, residual, xtx_inv, data[source])
    se_n = np.sqrt(np.diag(cov_naive))
    se_c = np.sqrt(np.clip(np.diag(cov_source), 0, None))
    j = list(xdf.columns).index(treatment)
    estimate = float(beta[j])
    se_naive = float(se_n[j])
    p_naive = float(2 * stats.t.sf(abs(estimate / se_naive), df=n - p))
    se_source = float(se_c[j])
    p_source = float(2 * stats.t.sf(abs(estimate / se_source), df=n_groups - 1))
    return AuditResult(
        outcome=outcome,
        estimate=estimate,
        se_naive=se_naive,
        p_naive=p_naive,
        se_source=se_source,
        p_source=p_source,
        se_inflation=se_source / se_naive,
        n_pairs=len(data),
        n_sources=data[source].nunique(),
    )


def wild_cluster_bootstrap(
    frame: pd.DataFrame,
    outcome: str,
    treatment: str = "treat",
    source: str = "source_id",
    covariates: Iterable[str] = (),
    reps: int = 999,
    seed: int = 20260711,
) -> dict:
    """Restricted Rademacher wild cluster bootstrap for one coefficient.

    The null model excludes the treatment coefficient.  One sign is drawn per
    source, preserving within-source dependence.  The returned p-value compares
    absolute source-clustered t statistics.
    """
    needed = [outcome, treatment, source, *covariates]
    data = frame.dropna(subset=needed).copy()
    xdf = _design(data, treatment, covariates)
    x0df = xdf.drop(columns=[treatment])
    x = xdf.to_numpy(float)
    x0 = x0df.to_numpy(float)
    y = data[outcome].astype(float).to_numpy()
    beta, residual_full, inv = _ols(y, x)
    cov, _ = _cluster_covariance(x, residual_full, inv, data[source])
    j = list(xdf.columns).index(treatment)
    t_obs = float(beta[j] / np.sqrt(cov[j, j]))
    beta0, residual, _ = _ols(y, x0)
    fitted = x0 @ beta0
    codes, groups = pd.factorize(data[source], sort=True)
    rng = np.random.default_rng(seed)
    exceed = 0
    valid = 0
    for _ in range(reps):
        weights = rng.choice((-1.0, 1.0), size=len(groups))
        y_star = fitted + residual * weights[codes]
        try:
            beta_star, resid_star, inv_star = _ols(y_star, x)
            cov_star, _ = _cluster_covariance(
                x, resid_star, inv_star, data[source]
            )
            t_star = float(beta_star[j] / np.sqrt(cov_star[j, j]))
        except (ValueError, np.linalg.LinAlgError):
            continue
        valid += 1
        exceed += abs(t_star) >= abs(t_obs)
    return {
        "outcome": outcome,
        "estimate": float(beta[j]),
        "t_observed": t_obs,
        "p_wild": (exceed + 1) / (valid + 1),
        "bootstrap_reps": valid,
        "n_sources": len(groups),
    }


def equal_source_test(
    frame: pd.DataFrame,
    outcome: str,
    treatment: str = "treat",
    source: str = "source_id",
    near_field: str | None = None,
) -> dict:
    """Give each source one weight by averaging its eligible recipients."""
    data = frame.copy()
    if near_field is not None:
        data = data.query(near_field)
    per_source = data.groupby([source, treatment], observed=True)[outcome].mean().reset_index()
    a = per_source.loc[per_source[treatment] == 1, outcome]
    b = per_source.loc[per_source[treatment] == 0, outcome]
    if len(a) < 2 or len(b) < 2:
        raise ValueError("each arm needs at least two sources")
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    welch = stats.ttest_ind(a, b, equal_var=False)
    return {
        "outcome": outcome,
        "difference": float(a.mean() - b.mean()),
        "p_mannwhitney": float(u.pvalue),
        "p_welch": float(welch.pvalue),
        "n_treated_sources": len(a),
        "n_control_sources": len(b),
    }


def leave_one_group_out(
    frame: pd.DataFrame,
    outcome: str,
    group: str = "mouse",
    **audit_kwargs,
) -> pd.DataFrame:
    """Repeat the source-aware audit after dropping each top-level group."""
    rows = []
    for value in sorted(frame[group].dropna().unique()):
        result = compare_replication_units(
            frame.loc[frame[group] != value], outcome=outcome, **audit_kwargs
        )
        rows.append({"dropped_group": value, **result.to_dict()})
    return pd.DataFrame(rows)


def benjamini_hochberg(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values in original order."""
    p = np.asarray(list(pvalues), dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out
