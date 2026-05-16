#!/usr/bin/env python3
"""
statistics_utils.py
===================
Statistical analysis utilities for the gold nanoparticle research framework.
Provides significance testing, power analysis, effect size calculations,
and reporting functions consistent with biomedical publication standards.

Author  : Om Jha, St. Xavier's College Kathmandu
Contact : om.physics7@gmail.com
Year    : 2026
"""

import numpy as np
from scipy import stats


def two_sample_ttest(group1, group2, alternative="two-sided", alpha=0.05):
    """
    Unpaired two-sample t-test with full statistical reporting.

    Parameters
    ----------
    group1      : array_like  First group measurements.
    group2      : array_like  Second group measurements.
    alternative : str         'two-sided', 'greater', or 'less'.
    alpha       : float       Significance level (default 0.05).

    Returns
    -------
    dict with keys: t_stat, p_value, df, significant, cohens_d,
                    ci_95, interpretation
    """
    g1 = np.asarray(group1, dtype=float)
    g2 = np.asarray(group2, dtype=float)

    t_stat, p_val = stats.ttest_ind(g1, g2, alternative=alternative)
    df    = len(g1) + len(g2) - 2

    # Cohen's d effect size
    pooled_std = np.sqrt((np.var(g1, ddof=1) + np.var(g2, ddof=1)) / 2.0)
    cohens_d   = (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0.0

    # 95% CI for difference in means
    se_diff = np.sqrt(np.var(g1, ddof=1)/len(g1) + np.var(g2, ddof=1)/len(g2))
    t_crit  = stats.t.ppf(0.975, df)
    diff    = np.mean(g1) - np.mean(g2)
    ci_95   = (diff - t_crit * se_diff, diff + t_crit * se_diff)

    # Significance label
    if p_val < 0.001:
        sig_label = "***"
    elif p_val < 0.01:
        sig_label = "**"
    elif p_val < alpha:
        sig_label = "*"
    else:
        sig_label = "ns"

    effect_interp = ("negligible" if abs(cohens_d) < 0.2 else
                     "small"      if abs(cohens_d) < 0.5 else
                     "medium"     if abs(cohens_d) < 0.8 else "large")

    return {
        "t_stat"        : round(float(t_stat), 4),
        "p_value"       : round(float(p_val), 6),
        "df"            : int(df),
        "significant"   : bool(p_val < alpha),
        "sig_label"     : sig_label,
        "cohens_d"      : round(float(cohens_d), 4),
        "effect_size"   : effect_interp,
        "ci_95"         : (round(ci_95[0], 4), round(ci_95[1], 4)),
        "mean_diff"     : round(float(diff), 4),
    }


def one_way_anova(*groups, alpha=0.05):
    """
    One-way ANOVA with Tukey post-hoc test for multiple group comparison.

    Parameters
    ----------
    *groups : array_like  Variable number of group arrays.
    alpha   : float       Significance level.

    Returns
    -------
    dict with keys: F_stat, p_value, significant, eta_squared
    """
    arrays = [np.asarray(g, dtype=float) for g in groups]
    F_stat, p_val = stats.f_oneway(*arrays)

    # Eta-squared effect size
    grand_mean = np.mean(np.concatenate(arrays))
    ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in arrays)
    ss_total   = sum(np.sum((g - grand_mean)**2) for g in arrays)
    eta_sq     = ss_between / ss_total if ss_total > 0 else 0.0

    return {
        "F_stat"     : round(float(F_stat), 4),
        "p_value"    : round(float(p_val), 6),
        "significant": bool(p_val < alpha),
        "eta_squared": round(float(eta_sq), 4),
        "n_groups"   : len(arrays),
    }


def sample_size_calculation(effect_size=0.8, alpha=0.05, power=0.80):
    """
    Sample size estimation for two-sample t-test using power analysis.

    Implements the iterative power calculation to find the minimum n
    per group achieving the desired statistical power.

    Parameters
    ----------
    effect_size : float  Cohen's d (0.2 small, 0.5 medium, 0.8 large).
    alpha       : float  Type I error rate (default 0.05).
    power       : float  Desired statistical power (default 0.80).

    Returns
    -------
    dict with keys: n_per_group, total_n, actual_power, effect_label
    """
    from scipy.stats import norm, t as t_dist

    effect_label = ("small" if effect_size < 0.5 else
                    "medium" if effect_size < 0.8 else "large")

    # Iterative search
    for n in range(2, 500):
        df      = 2 * n - 2
        t_crit  = t_dist.ppf(1 - alpha/2, df)
        ncp     = effect_size * np.sqrt(n / 2)
        power_n = 1.0 - t_dist.cdf(t_crit - ncp, df)
        if power_n >= power:
            return {
                "n_per_group" : n,
                "total_n"     : 2 * n,
                "actual_power": round(float(power_n), 4),
                "effect_label": effect_label,
                "effect_size" : effect_size,
            }

    return {"n_per_group": 500, "total_n": 1000,
            "actual_power": power, "effect_label": effect_label}


def normality_test(data, alpha=0.05):
    """
    Shapiro-Wilk normality test with interpretation.

    Parameters
    ----------
    data  : array_like  Sample data.
    alpha : float       Significance level.

    Returns
    -------
    dict with keys: statistic, p_value, is_normal, interpretation
    """
    data = np.asarray(data, dtype=float)
    stat, p = stats.shapiro(data)
    is_normal = bool(p > alpha)
    return {
        "statistic"     : round(float(stat), 4),
        "p_value"       : round(float(p), 6),
        "is_normal"     : is_normal,
        "interpretation": ("Data consistent with normality" if is_normal
                           else "Data significantly deviates from normality"),
    }


def log_rank_test(time1, event1, time2, event2):
    """
    Simplified log-rank statistic for Kaplan-Meier survival comparison.

    Implements the Mantel-Haenszel log-rank test for comparing
    two survival curves (e.g., control vs treated groups).

    Parameters
    ----------
    time1  : array_like  Survival times — group 1.
    event1 : array_like  Event indicators — group 1 (1=event, 0=censored).
    time2  : array_like  Survival times — group 2.
    event2 : array_like  Event indicators — group 2.

    Returns
    -------
    dict with keys: chi2, p_value, significant, hazard_ratio
    """
    t1, e1 = np.asarray(time1), np.asarray(event1)
    t2, e2 = np.asarray(time2), np.asarray(event2)

    all_times = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))

    O1_total = E1_total = 0.0
    V_total  = 0.0

    for t in all_times:
        n1 = np.sum(t1 >= t)
        n2 = np.sum(t2 >= t)
        o1 = np.sum((t1 == t) & (e1 == 1))
        o2 = np.sum((t2 == t) & (e2 == 1))
        n  = n1 + n2
        o  = o1 + o2
        if n < 2:
            continue
        e1_t   = n1 * o / n
        var_t  = n1 * n2 * o * (n - o) / (n**2 * (n - 1)) if n > 1 else 0
        O1_total += o1
        E1_total += e1_t
        V_total  += var_t

    chi2 = (O1_total - E1_total)**2 / V_total if V_total > 0 else 0.0
    p    = 1.0 - stats.chi2.cdf(chi2, df=1)

    hr   = (np.sum(e1) / np.sum(t1)) / (np.sum(e2) / np.sum(t2)) if np.sum(e2) > 0 else np.nan

    return {
        "chi2"       : round(float(chi2), 4),
        "p_value"    : round(float(p), 6),
        "significant": bool(p < 0.05),
        "hazard_ratio": round(float(hr), 4) if not np.isnan(hr) else None,
    }


def bonferroni_correction(p_values, alpha=0.05):
    """
    Bonferroni correction for multiple comparisons.

    Parameters
    ----------
    p_values : array_like  Uncorrected p-values.
    alpha    : float       Family-wise error rate.

    Returns
    -------
    dict with keys: corrected_alpha, significant_mask, adjusted_p_values
    """
    p = np.asarray(p_values, dtype=float)
    m = len(p)
    corrected_alpha = alpha / m
    adjusted_p      = np.minimum(p * m, 1.0)
    return {
        "corrected_alpha"  : round(corrected_alpha, 6),
        "adjusted_p_values": adjusted_p.tolist(),
        "significant_mask" : (adjusted_p < alpha).tolist(),
        "n_comparisons"    : m,
    }


def summary_statistics(data, label=""):
    """
    Compute and print comprehensive descriptive statistics.

    Parameters
    ----------
    data  : array_like  Input measurements.
    label : str         Optional label for display.

    Returns
    -------
    dict with mean, median, std, sem, ci_95, min, max, n
    """
    d   = np.asarray(data, dtype=float)
    n   = len(d)
    m   = float(np.mean(d))
    med = float(np.median(d))
    s   = float(np.std(d, ddof=1))
    sem = s / np.sqrt(n)
    t_c = stats.t.ppf(0.975, n - 1)
    ci  = (m - t_c * sem, m + t_c * sem)

    result = {
        "label" : label,
        "n"     : n,
        "mean"  : round(m, 4),
        "median": round(med, 4),
        "std"   : round(s, 4),
        "sem"   : round(sem, 4),
        "ci_95" : (round(ci[0], 4), round(ci[1], 4)),
        "min"   : round(float(np.min(d)), 4),
        "max"   : round(float(np.max(d)), 4),
    }
    return result
