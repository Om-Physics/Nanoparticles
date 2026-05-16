#!/usr/bin/env python3
"""
quick_start.py
==============
Minimal working example demonstrating the core computational analyses.

Run from repository root:
    python examples/quick_start.py

Author  : Om Jha, St. Xavier's College Kathmandu
Contact : om.physics7@gmail.com
Year    : 2026
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

# ── imports ───────────────────────────────────────────────────────────────────
from src.models import (
    two_compartment_pk,
    epr_tumor_accumulation,
    calculate_pk_parameters,
    photothermal_heating,
    cem43_dose,
    ph_dependent_release,
    korsmeyer_peppas_release,
    fit_ic50,
    four_parameter_logistic,
)
from src.statistics_utils import (
    two_sample_ttest,
    one_way_anova,
    sample_size_calculation,
    normality_test,
    summary_statistics,
    bonferroni_correction,
)
from src.generate_figures import (
    figure_biodistribution,
    figure_toxicity,
    figure_mechanistic_controls,
    figure_drug_loading,
    figure_photothermal,
    figure_pharmacokinetics,
)

DIVIDER = "=" * 68


def section(title):
    print(f"\n{DIVIDER}")
    print(f"  {title}")
    print(DIVIDER)


# ══════════════════════════════════════════════════════════════════════════════
# 1. PHARMACOKINETIC PROFILES
# ══════════════════════════════════════════════════════════════════════════════

section("1. Pharmacokinetic Profiles")

t = np.linspace(0, 168, 2000)   # 0 to 7 days, 2000 points

# Free doxorubicin — rapid mono-exponential clearance
C_free    = 100 * np.exp(-0.35 * t)

# Passive PEG-AuNP-DOX — bi-exponential, t½β = 46 h
C_passive = two_compartment_pk(t, A=65, alpha=0.277, B=35, beta=0.015)

# Active C225-AuNP-DOX — bi-exponential, t½β = 58 h
C_active  = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)

pk_free    = calculate_pk_parameters(t, C_free)
pk_passive = calculate_pk_parameters(t, C_passive)
pk_active  = calculate_pk_parameters(t, C_active)

for label, pk in [("Free DOX", pk_free),
                  ("Passive AuNP-DOX", pk_passive),
                  ("Active AuNP-DOX",  pk_active)]:
    print(f"\n{label}")
    print(f"  Cmax     = {pk['Cmax']:.2f} %ID/mL")
    print(f"  AUC      = {pk['AUC']:.1f} %ID·h/mL")
    print(f"  t½ (term)= {pk['t_half']:.1f} h")

fold_auc = pk_active["AUC"] / pk_free["AUC"]
print(f"\nAUC fold-improvement (active vs free): {fold_auc:.1f}×")


# ══════════════════════════════════════════════════════════════════════════════
# 2. TUMOUR ACCUMULATION — EPR MODEL
# ══════════════════════════════════════════════════════════════════════════════

section("2. EPR-Mediated Tumour Accumulation")

plasma_passive = (65, 0.277, 35, 0.015)
plasma_active  = (70, 0.277, 30, 0.012)

C_tum_passive = epr_tumor_accumulation(
    t, k_in=0.085, k_out=0.015, plasma_params=plasma_passive
)
C_tum_active  = epr_tumor_accumulation(
    t, k_in=0.110, k_out=0.012, plasma_params=plasma_active
)

peak_p_val = C_tum_passive.max()
peak_p_t   = t[np.argmax(C_tum_passive)]
peak_a_val = C_tum_active.max()
peak_a_t   = t[np.argmax(C_tum_active)]

print(f"\nPassive AuNP-DOX  peak = {peak_p_val:.2f} %ID/g  at {peak_p_t:.1f} h")
print(f"Active  AuNP-DOX  peak = {peak_a_val:.2f} %ID/g  at {peak_a_t:.1f} h")
print(f"Enhancement             = {peak_a_val/peak_p_val*100 - 100:.1f}%")


# ══════════════════════════════════════════════════════════════════════════════
# 3. PH-RESPONSIVE DRUG RELEASE
# ══════════════════════════════════════════════════════════════════════════════

section("3. pH-Responsive Drug Release Kinetics")

t_rel = np.linspace(0, 72, 720)

for ph, label in [(7.4, "Blood pH 7.4"),
                  (6.5, "Tumour pH 6.5"),
                  (5.0, "Endosome pH 5.0")]:
    released = ph_dependent_release(t_rel, pH=ph)
    print(f"  {label:20s}  at 72 h: {released[-1]:.1f}%")

# Korsmeyer-Peppas fit example
kp_release = korsmeyer_peppas_release(t_rel, k=0.018, n=0.52)
print(f"\nKorsmeyer-Peppas (pH 7.4) release exponent n = 0.52 (anomalous transport)")
print(f"Release at 24 h: {kp_release[np.searchsorted(t_rel, 24)]:.1f}%")


# ══════════════════════════════════════════════════════════════════════════════
# 4. PHOTOTHERMAL THERAPY
# ══════════════════════════════════════════════════════════════════════════════

section("4. Photothermal Therapy — CEM43 Thermal Dose")

t_ptt  = np.linspace(0, 600, 600)    # 10 minutes in seconds
T_ptt  = photothermal_heating(t_ptt, eta=0.99, I=1.5, tau=90, T_body=37.0)
cem    = cem43_dose(T_ptt, dt=1.0)

print(f"\nLaser: 808 nm, 1.5 W/cm², 10 min exposure")
print(f"Max temperature:         {T_ptt.max():.1f} °C")
print(f"Steady-state ΔT:         {T_ptt.max() - 37:.1f} °C above baseline")
print(f"CEM43 at end of exposure: {cem[-1]:.1f} min-equivalent")
print(f"Necrosis threshold (60): {'EXCEEDED' if cem[-1] > 60 else 'Not reached'}")


# ══════════════════════════════════════════════════════════════════════════════
# 5. IC50 DETERMINATION
# ══════════════════════════════════════════════════════════════════════════════

section("5. IC50 Dose-Response Analysis (A549 Cells)")

doses_log = np.logspace(0, 3, 20)   # 1–1000 µg/mL

# Simulated viability data with noise
np.random.seed(42)
viab_free = (four_parameter_logistic(doses_log, 0, 100, 52, -1.8)
             + np.random.normal(0, 2, 20))
viab_np   = (four_parameter_logistic(doses_log, 0, 100, 28, -2.1)
             + np.random.normal(0, 2, 20))

ic50_free, _, _, r2_free = fit_ic50(doses_log, viab_free)
ic50_np,   _, _, r2_np   = fit_ic50(doses_log, viab_np)

print(f"\nFree doxorubicin:   IC50 = {ic50_free:.1f} µg/mL  (R² = {r2_free:.3f})")
print(f"AuNP-DOX (active):  IC50 = {ic50_np:.1f} µg/mL  (R² = {r2_np:.3f})")
print(f"Fold improvement:   {ic50_free/ic50_np:.2f}×")


# ══════════════════════════════════════════════════════════════════════════════
# 6. STATISTICAL ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

section("6. Statistical Analysis")

# Simulated tumour accumulation replicates (n=8 per group)
np.random.seed(0)
passive_rep = np.random.normal(4.8, 0.35, 8)
active_rep  = np.random.normal(8.5, 0.40, 8)

# Normality check
norm_p = normality_test(passive_rep)
norm_a = normality_test(active_rep)
print(f"\nShapiro-Wilk normality: passive p={norm_p['p_value']}, "
      f"active p={norm_a['p_value']}")

# t-test
ttest = two_sample_ttest(passive_rep, active_rep)
print(f"\nTumour accumulation comparison (passive vs active):")
print(f"  t({ttest['df']}) = {ttest['t_stat']},  p = {ttest['p_value']}  "
      f"{ttest['sig_label']}")
print(f"  Cohen's d = {ttest['cohens_d']}  ({ttest['effect_size']} effect)")
print(f"  95% CI for difference: {ttest['ci_95']}")

# One-way ANOVA across all groups
control  = np.random.normal(1.2, 0.2, 8)
free_dox = np.random.normal(3.1, 0.4, 8)
anova    = one_way_anova(control, free_dox, passive_rep, active_rep)
print(f"\nOne-way ANOVA (4 groups):")
print(f"  F = {anova['F_stat']},  p = {anova['p_value']},  η² = {anova['eta_squared']}")

# Bonferroni correction for 6 pairwise comparisons
raw_p      = [0.041, 0.0003, 0.0001, 0.028, 0.012, 0.0008]
corrected  = bonferroni_correction(raw_p)
print(f"\nBonferroni correction (6 comparisons):")
print(f"  Corrected α = {corrected['corrected_alpha']}")
print(f"  Significant: {corrected['significant_mask']}")

# Power analysis
power = sample_size_calculation(effect_size=0.8, alpha=0.05, power=0.80)
print(f"\nPower analysis (Cohen's d=0.8, α=0.05, power=0.80):")
print(f"  Required n per group: {power['n_per_group']}")
print(f"  Total animals:        {power['total_n']}")
print(f"  Achieved power:       {power['actual_power']:.3f}")

# Summary statistics
stats_a = summary_statistics(active_rep, label="Active AuNP-DOX tumour accumulation")
print(f"\n{stats_a['label']}:")
print(f"  Mean ± SD: {stats_a['mean']} ± {stats_a['std']} %ID/g")
print(f"  95% CI:    {stats_a['ci_95']}")


# ══════════════════════════════════════════════════════════════════════════════
# 7. GENERATE ALL FIGURES
# ══════════════════════════════════════════════════════════════════════════════

section("7. Generating All Publication Figures")

figure_functions = [
    ("Figure 1 — Biodistribution",         figure_biodistribution),
    ("Figure 2 — Toxicity Panel",          figure_toxicity),
    ("Figure 3 — Mechanistic Controls",    figure_mechanistic_controls),
    ("Figure 4 — Drug Loading",            figure_drug_loading),
    ("Figure 5 — Photothermal Therapy",    figure_photothermal),
    ("Figure 6 — Pharmacokinetics",        figure_pharmacokinetics),
]

for name, fn in figure_functions:
    print(f"\n  Generating {name}...")
    fn()

print(f"\n{DIVIDER}")
print("  All figures saved to ./figures/ at 300 DPI")
print(f"{DIVIDER}\n")
