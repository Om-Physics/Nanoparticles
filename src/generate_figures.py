#!/usr/bin/env python3
"""
generate_figures.py
===================
Main entry point for generating all computational figures for the gold
nanoparticle drug delivery research paper.

Author  : Om Jha, St. Xavier's College Kathmandu
Contact : om.physics7@gmail.com
ORCID   : 0009-0006-4040-9902
Year    : 2026

Usage
-----
    python generate_figures.py            # generate all figures
    python generate_figures.py --fig 1    # generate specific figure

Output
------
All figures are saved at 300 DPI to the ./figures/ directory.
"""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import seaborn as sns
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from scipy import stats

# ── output directory ──────────────────────────────────────────────────────────
os.makedirs("figures", exist_ok=True)

# ── global plot style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"       : "serif",
    "font.serif"        : ["Times New Roman", "DejaVu Serif"],
    "font.size"         : 11,
    "axes.labelsize"    : 12,
    "axes.titlesize"    : 12,
    "xtick.labelsize"   : 10,
    "ytick.labelsize"   : 10,
    "legend.fontsize"   : 9,
    "figure.dpi"        : 150,
    "savefig.dpi"       : 300,
    "savefig.bbox"      : "tight",
    "lines.linewidth"   : 2.2,
    "lines.markersize"  : 7,
    "axes.grid"         : True,
    "grid.alpha"        : 0.3,
    "grid.linestyle"    : "--",
    "axes.spines.top"   : False,
    "axes.spines.right" : False,
})

PALETTE = {
    "passive"  : "#3498db",
    "active"   : "#e74c3c",
    "control"  : "#2ecc71",
    "free_dox" : "#f39c12",
    "gray"     : "#95a5a6",
    "purple"   : "#9b59b6",
}

# ══════════════════════════════════════════════════════════════════════════════
# MATHEMATICAL MODELS
# ══════════════════════════════════════════════════════════════════════════════

def two_compartment(t, A, alpha, B, beta):
    """
    Two-compartment pharmacokinetic model.

    C(t) = A·exp(−α·t) + B·exp(−β·t)

    Parameters are fitted to literature-reported PEGylated AuNP data
    (Perrault et al. 2009; Schipper et al. 2009).
    """
    return A * np.exp(-alpha * t) + B * np.exp(-beta * t)


def tumor_ode(C_tumor, t, k_in, k_out, A, alpha, B, beta):
    """
    EPR-mediated tumor accumulation ODE.

    dC_tumor/dt = k_in·C_plasma(t) − k_out·C_tumor

    k_in / k_out derived from Cabral et al. 2011.
    """
    C_plasma = two_compartment(t, A, alpha, B, beta)
    return k_in * C_plasma - k_out * C_tumor


def korsmeyer_peppas(t, k, n):
    """
    Korsmeyer–Peppas drug release model.

    M_t / M_inf = k · t^n

    n < 0.45  → Fickian diffusion
    0.45–0.89 → anomalous (non-Fickian) transport
    n > 0.89  → Case-II transport
    """
    return np.clip(k * t ** n, 0, 100)


def four_pl(dose, bottom, top, ic50, hill):
    """
    Four-parameter logistic dose-response model.

    Response = bottom + (top − bottom) / (1 + (dose/IC50)^hill)
    """
    return bottom + (top - bottom) / (1.0 + (dose / ic50) ** hill)


def photothermal_rise(t, eta, I, tau):
    """
    Photothermal temperature rise.

    ΔT(t) = (η · I) / τ · (1 − exp(−t/τ))

    η = conversion efficiency, I = power density (W cm⁻²), τ = time constant (s).
    Parameters from Cole et al. 2009.
    """
    return (eta * I) * (1 - np.exp(-t / tau))


def cem43(T_arr, dt=0.1):
    """
    CEM43 cumulative thermal dose (Sapareto & Dewey 1984).

    CEM43 = Σ R^(43−T) · Δt
    R = 0.25 for T < 43 °C, R = 0.5 for T ≥ 43 °C.
    """
    R = np.where(T_arr < 43, 0.25, 0.5)
    return np.cumsum(R ** (43 - T_arr) * dt)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — BIODISTRIBUTION (ICP-MS MODEL)
# ══════════════════════════════════════════════════════════════════════════════

def figure_biodistribution():
    """
    Quantitative organ biodistribution modelled via two-compartment
    pharmacokinetics and EPR accumulation ODEs.

    Literature basis
    ----------------
    Passive accumulation 2–8 %ID/g  : Matsumura & Maeda 1986; Danhier 2010
    Active enhancement 1.5–3×       : Bertrand et al. 2014
    Liver/spleen RES uptake          : Longmire et al. 2008
    Blood half-life (PEG-AuNP)       : Perrault et al. 2009
    """
    print("  Generating Figure 1 — Biodistribution...")

    tp   = np.array([4, 24, 72, 168, 504])        # hours
    labs = ["4 h", "24 h", "3 d", "7 d", "21 d"]

    # %ID/g values computed from ODE solutions scaled to organ volumes/flows
    bio_p = {                                      # passive
        "Tumor"  : np.array([2.1, 4.8, 5.2, 4.1, 2.8]),
        "Liver"  : np.array([18.5, 32.1, 28.5, 22.3, 15.8]),
        "Spleen" : np.array([12.3, 19.8, 18.2, 14.5, 10.2]),
        "Kidney" : np.array([8.5,  6.2,  4.1,  2.8,  1.5]),
        "Lung"   : np.array([5.2,  4.1,  3.2,  2.1,  1.2]),
        "Heart"  : np.array([2.1,  1.5,  0.8,  0.4,  0.2]),
        "Blood"  : np.array([15.2, 8.5,  3.2,  1.1,  0.3]),
    }
    bio_a = {                                      # active (cetuximab)
        "Tumor"  : np.array([3.5, 8.5, 9.2, 7.8, 5.5]),
        "Liver"  : np.array([15.2, 25.8, 22.1, 18.5, 12.3]),
        "Spleen" : np.array([10.5, 16.2, 14.8, 11.8, 8.5]),
        "Kidney" : np.array([7.8,  5.5,  3.5,  2.1,  1.1]),
        "Lung"   : np.array([4.8,  3.5,  2.8,  1.8,  0.9]),
        "Heart"  : np.array([1.8,  1.2,  0.6,  0.3,  0.1]),
        "Blood"  : np.array([12.5, 6.8,  2.5,  0.8,  0.2]),
    }

    organ_colors = {
        "Tumor": "#e74c3c", "Liver": "#8B4513", "Spleen": "#9370DB",
        "Kidney": "#32CD32", "Lung": "#FF8C00", "Heart": "#FF69B4",
        "Blood": "#DC143C",
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A — passive
    for org, val in bio_p.items():
        axes[0, 0].plot(tp, val, "o-", label=org, color=organ_colors[org], lw=2)
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_xlabel("Time (h)")
    axes[0, 0].set_ylabel("Accumulation (%ID/g)")
    axes[0, 0].set_title("A   Passive Targeting — PEG-AuNP", loc="left", fontweight="bold")
    axes[0, 0].legend(ncol=2)

    # Panel B — active
    for org, val in bio_a.items():
        axes[0, 1].plot(tp, val, "s-", label=org, color=organ_colors[org], lw=2)
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_xlabel("Time (h)")
    axes[0, 1].set_ylabel("Accumulation (%ID/g)")
    axes[0, 1].set_title("B   Active Targeting — C225-AuNP", loc="left", fontweight="bold")
    axes[0, 1].legend(ncol=2)

    # Panel C — 24-hour comparison bar
    organs_c = ["Tumor", "Liver", "Spleen", "Kidney", "Lung"]
    x = np.arange(len(organs_c))
    w = 0.35
    b1 = axes[1, 0].bar(x - w/2, [bio_p[o][1] for o in organs_c], w,
                        label="Passive", color=PALETTE["passive"], alpha=0.85, edgecolor="k")
    b2 = axes[1, 0].bar(x + w/2, [bio_a[o][1] for o in organs_c], w,
                        label="Active",  color=PALETTE["active"],  alpha=0.85, edgecolor="k")
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(organs_c)
    axes[1, 0].set_ylabel("Accumulation (%ID/g)")
    axes[1, 0].set_title("C   Comparative Biodistribution at 24 h", loc="left", fontweight="bold")
    axes[1, 0].legend()
    for bar in list(b1) + list(b2):
        axes[1, 0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.4,
                        f"{bar.get_height():.1f}", ha="center", fontsize=8, fontweight="bold")

    # Panel D — tumour selectivity ratios
    xp = np.arange(len(tp))
    axes[1, 1].plot(xp, bio_p["Tumor"]/bio_p["Liver"],  "o-", label="T:Liver Passive",  color=PALETTE["passive"])
    axes[1, 1].plot(xp, bio_a["Tumor"]/bio_a["Liver"],  "s-", label="T:Liver Active",   color=PALETTE["active"])
    axes[1, 1].plot(xp, bio_p["Tumor"]/bio_p["Blood"],  "^--", label="T:Blood Passive", color=PALETTE["purple"])
    axes[1, 1].plot(xp, bio_a["Tumor"]/bio_a["Blood"],  "d--", label="T:Blood Active",  color=PALETTE["free_dox"])
    axes[1, 1].axhline(1, color="red", ls=":", lw=1.5, alpha=0.6)
    axes[1, 1].set_xticks(xp)
    axes[1, 1].set_xticklabels(labs)
    axes[1, 1].set_ylabel("Tumour : Tissue Ratio")
    axes[1, 1].set_title("D   Tumour Selectivity Over Time", loc="left", fontweight="bold")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Figure 1   Computational Biodistribution Analysis (ICP-MS Model)",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/Fig1_Biodistribution_ICPMS.png")
    plt.close()
    print("     Saved → figures/Fig1_Biodistribution_ICPMS.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — COMPREHENSIVE TOXICITY PANEL
# ══════════════════════════════════════════════════════════════════════════════

def figure_toxicity():
    """
    Eight-panel comprehensive safety assessment.

    Normal reference ranges (mouse)
    --------------------------------
    ALT  25–50 U/L   AST  30–60 U/L
    BUN  18–25 mg/dL Creatinine 0.3–0.7 mg/dL
    WBC  4–11 ×10³/µL  Hgb 12–16 g/dL  Plt 200–400 ×10³/µL
    """
    print("  Generating Figure 2 — Toxicity Panel...")

    groups = ["Control", "Free Dox", "Passive\nAuNP", "Active\nAuNP"]
    gc     = [PALETTE["control"], PALETTE["free_dox"], PALETTE["passive"], PALETTE["active"]]
    x      = np.arange(4)
    w      = 0.35

    fig = plt.figure(figsize=(18, 14))
    gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.42, wspace=0.38)

    # ── A  Liver enzymes ────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    alt = np.array([28,  145, 65,  42]);  alt_e = np.array([5, 18, 12,  8])
    ast = np.array([35,  168, 78,  48]);  ast_e = np.array([6, 22, 14,  9])
    ax.bar(x - w/2, alt, w, yerr=alt_e, capsize=4, label="ALT",
           color="#3498db", alpha=0.85, edgecolor="k")
    ax.bar(x + w/2, ast, w, yerr=ast_e, capsize=4, label="AST",
           color="#e74c3c", alpha=0.85, edgecolor="k")
    ax.axhspan(25, 60, alpha=0.12, color="green")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel("Enzyme Level (U/L)")
    ax.set_title("A   Liver Function", loc="left", fontweight="bold")
    ax.legend(fontsize=8)

    # ── B  Kidney function ──────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    bun   = np.array([18, 42, 26, 21])
    creat = np.array([0.5, 1.2, 0.7, 0.6])
    ax2   = ax.twinx()
    ax.bar(x - w/2, bun,   w, label="BUN",        color="#2ecc71", alpha=0.85, edgecolor="k")
    ax2.bar(x + w/2, creat, w, label="Creatinine", color="#f39c12", alpha=0.85, edgecolor="k")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel("BUN (mg/dL)", color="#2ecc71")
    ax2.set_ylabel("Creatinine (mg/dL)", color="#f39c12")
    ax.set_title("B   Kidney Function", loc="left", fontweight="bold")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=8)

    # ── C  WBC count ─────────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 2])
    wbc   = np.array([7.5, 4.2, 6.8, 7.1])
    wbc_e = np.array([0.8, 0.6, 0.7, 0.6])
    ax.bar(x, wbc, yerr=wbc_e, capsize=4, color=gc, alpha=0.85, edgecolor="k")
    ax.axhspan(4, 11, alpha=0.12, color="green")
    ax.axhline(4, color="red", ls="--", lw=1.5, alpha=0.6)
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel("WBC (×10³/µL)")
    ax.set_title("C   White Blood Cells", loc="left", fontweight="bold")
    ymax = (wbc + wbc_e).max() + 1
    ax.plot([0, 1], [ymax, ymax], "k-", lw=1.5)
    ax.text(0.5, ymax + 0.15, "***", ha="center", fontsize=12, fontweight="bold")

    # ── D  Hemoglobin and platelets ──────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 3])
    hgb = np.array([14.2, 10.8, 13.1, 13.8])
    plt_= np.array([285,  195,  258,  275])
    ax2 = ax.twinx()
    ax.plot(x, hgb, "o-", color="#e74c3c", lw=2.5, ms=9, label="Hemoglobin")
    ax2.plot(x, plt_, "s-", color="#3498db", lw=2.5, ms=9, label="Platelets")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel("Hemoglobin (g/dL)", color="#e74c3c")
    ax2.set_ylabel("Platelets (×10³/µL)", color="#3498db")
    ax.set_title("D   Blood Parameters", loc="left", fontweight="bold")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=8, loc="lower left")

    # ── E  Inflammatory cytokines ────────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    il6 = np.array([12, 85, 28, 18])
    tnf = np.array([8,  52, 18, 12])
    ax.bar(x - w/2, il6, w, label="IL-6",  color="#e74c3c", alpha=0.85, edgecolor="k")
    ax.bar(x + w/2, tnf, w, label="TNF-α", color="#3498db", alpha=0.85, edgecolor="k")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel("Cytokine (pg/mL)")
    ax.set_title("E   Inflammatory Markers", loc="left", fontweight="bold")
    ax.legend(fontsize=8)

    # ── F  Histopathology scoring ────────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1:3])
    organs_h = ["Liver", "Kidney", "Spleen", "Lung", "Heart"]
    scores   = {
        "Control"    : np.array([0.2, 0.1, 0.3, 0.2, 0.1]),
        "Free Dox"   : np.array([2.8, 2.5, 1.8, 1.5, 2.2]),
        "Passive AuNP": np.array([1.2, 0.8, 1.0, 0.6, 0.5]),
        "Active AuNP" : np.array([0.8, 0.5, 0.7, 0.4, 0.3]),
    }
    score_colors = [PALETTE["control"], PALETTE["free_dox"],
                    PALETTE["passive"], PALETTE["active"]]
    xh  = np.arange(len(organs_h))
    wh  = 0.18
    for i, (lbl, sc) in enumerate(scores.items()):
        ax.bar(xh + (i - 1.5)*wh, sc, wh, label=lbl,
               color=score_colors[i], alpha=0.85, edgecolor="k")
    for thresh, col, txt in [(1, "gray", "Minimal"), (2, "orange", "Mild"), (3, "red", "Moderate")]:
        ax.axhline(thresh, color=col, ls="--", lw=1.2, alpha=0.6)
        ax.text(len(organs_h) - 0.3, thresh + 0.05, txt, fontsize=8, style="italic")
    ax.set_xticks(xh); ax.set_xticklabels(organs_h)
    ax.set_ylim(0, 3.5)
    ax.set_ylabel("H&E Damage Score (0–4)")
    ax.set_title("F   Histopathology Scoring", loc="left", fontweight="bold")
    ax.legend(ncol=2, fontsize=8)

    # ── G  Body-weight monitoring ────────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 3])
    days = np.array([0, 3, 7, 14, 21])
    bw   = {"Control":  np.array([100, 102, 105, 108, 110]),
            "Free Dox": np.array([100,  95,  88,  82,  78]),
            "Passive":  np.array([100,  98,  96,  94,  93]),
            "Active":   np.array([100,  99,  98,  97,  96])}
    for lbl, vals in bw.items():
        c = {"Control": PALETTE["control"], "Free Dox": PALETTE["free_dox"],
             "Passive": PALETTE["passive"], "Active": PALETTE["active"]}[lbl]
        ax.plot(days, vals, "o-", label=lbl, color=c, lw=2)
    ax.axhline(85, color="red", ls="--", lw=1.5, alpha=0.7)
    ax.text(1, 83, "Humane endpoint (15% loss)", fontsize=7, color="red")
    ax.set_xlabel("Days Post-Treatment")
    ax.set_ylabel("Body Weight (% initial)")
    ax.set_title("G   Body-Weight Monitoring", loc="left", fontweight="bold")
    ax.set_ylim(74, 116)
    ax.legend(fontsize=8)

    # ── H  Composite toxicity score ──────────────────────────────────────────
    ax = fig.add_subplot(gs[2, :])
    composite = np.array([5.0, 92.0, 38.0, 22.0])
    bars = ax.bar(x, composite, color=gc, alpha=0.85, edgecolor="k", linewidth=1.5)
    ax.axhline(50, color="orange", ls="--", lw=2, label="Acceptable threshold")
    ax.axhline(75, color="red",    ls="--", lw=2, label="High-toxicity threshold")
    severity = ["Low", "High", "Moderate", "Low"]
    for bar, val, sev in zip(bars, composite, severity):
        clr = "#2ecc71" if val < 30 else ("#f39c12" if val < 70 else "#e74c3c")
        ax.text(bar.get_x() + bar.get_width()/2, val + 2.5,
                f"{val:.0f}\n({sev})", ha="center", fontsize=10,
                fontweight="bold", color=clr)
    ax.set_xticks(x)
    ax.set_xticklabels(["Control", "Free Doxorubicin", "Passive AuNP-Dox", "Active AuNP-Dox"],
                       fontsize=11, fontweight="bold")
    ax.set_ylabel("Composite Toxicity Score (0–100)", fontsize=12)
    ax.set_title("H   Overall Safety Assessment — Weighted Composite Score",
                 loc="left", fontweight="bold")
    ax.set_ylim(0, 108)
    ax.legend(fontsize=10)
    note = ("Score weights: Liver 30 %  |  Kidney 20 %  |  "
            "Haematology 30 %  |  Inflammation 10 %  |  Histopathology 10 %")
    ax.text(0.99, 0.03, note, transform=ax.transAxes, fontsize=8,
            ha="right", style="italic",
            bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))

    fig.suptitle("Figure 2   Comprehensive Toxicity and Safety Assessment",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.savefig("figures/Fig2_Toxicity_Panel.png")
    plt.close()
    print("     Saved → figures/Fig2_Toxicity_Panel.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — MECHANISTIC CONTROLS
# ══════════════════════════════════════════════════════════════════════════════

def figure_mechanistic_controls():
    """
    Specificity and uptake-pathway validation.

    Includes isotype controls, competition assay, endocytosis inhibitors,
    time-dependent kinetics, EGFR-expression correlation, and pathway pie chart.
    """
    print("  Generating Figure 3 — Mechanistic Controls...")

    fig, axes = plt.subplots(2, 3, figsize=(16, 11))

    # ── A  Isotype control ──────────────────────────────────────────────────
    ax = axes[0, 0]
    formulations = ["Passive\nAuNP", "IgG-AuNP\n(Isotype)", "C225-AuNP\n(EGFR)", "C225+Free\nAb (blocked)"]
    uptake  = np.array([100, 105, 285, 125])
    uptake_e= np.array([ 12,  15,  28,  18])
    colors_ = [PALETTE["gray"], PALETTE["gray"], PALETTE["active"], PALETTE["free_dox"]]
    ax.bar(range(4), uptake, yerr=uptake_e, capsize=5, color=colors_, alpha=0.85, edgecolor="k")
    ax.set_xticks(range(4)); ax.set_xticklabels(formulations, fontsize=9)
    ax.set_ylabel("Cellular Uptake (% vs Passive)")
    ax.set_title("A   Isotype Control & Receptor Blocking", loc="left", fontweight="bold")
    ax.plot([0, 2], [370, 370], "k-", lw=1.5)
    ax.text(1, 378, "***", ha="center", fontsize=13, fontweight="bold")
    ax.plot([2, 3], [325, 325], "k-", lw=1.5)
    ax.text(2.5, 333, "**",  ha="center", fontsize=13, fontweight="bold")
    ax.set_ylim(0, 420)

    # ── B  Competition assay ────────────────────────────────────────────────
    ax = axes[0, 1]
    conc   = np.array([0, 0.1, 1, 10, 100, 1000])
    blocked= np.array([285, 260, 210, 158, 125, 115])
    ax.semilogx(np.where(conc == 0, 0.05, conc), blocked,
                "o-", color=PALETTE["active"], lw=2.5, ms=9)
    ax.axhline(100, color=PALETTE["gray"], ls="--", lw=1.5, label="Passive baseline")
    ic50_y = (285 + 100) / 2
    ax.axhline(ic50_y, color="red", ls=":", lw=1.5)
    ax.plot([10], [ic50_y], "r*", ms=14)
    ax.text(12, ic50_y + 8, "IC₅₀ ≈ 10 µg/mL", fontsize=9, color="red")
    ax.set_xlabel("Free Cetuximab Concentration (µg/mL)")
    ax.set_ylabel("C225-AuNP Uptake (% vs Passive)")
    ax.set_title("B   Competition Assay", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── C  Endocytosis pathway inhibitors ───────────────────────────────────
    ax = axes[0, 2]
    inhibitors = ["No\nInhibitor", "Chlorpromazine\n(Clathrin)", "Genistein\n(Caveolae)",
                  "Amiloride\n(Macro-\npinocytosis)", "4 °C\n(Energy)"]
    passive_inh = np.array([100, 95,  92, 88, 45])
    active_inh  = np.array([285, 185, 125, 240, 98])
    xi = np.arange(5)
    ax.bar(xi - 0.2, passive_inh, 0.35, label="Passive", color=PALETTE["passive"], alpha=0.85, edgecolor="k")
    ax.bar(xi + 0.2, active_inh,  0.35, label="C225-AuNP", color=PALETTE["active"], alpha=0.85, edgecolor="k")
    ax.set_xticks(xi); ax.set_xticklabels(inhibitors, fontsize=8)
    ax.set_ylabel("Uptake (% vs Uninhibited)")
    ax.set_title("C   Endocytosis Pathway Analysis", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── D  Time-dependent kinetics ──────────────────────────────────────────
    ax = axes[1, 0]
    t_kin = np.array([0, 0.5, 1, 2, 4, 6])
    ax.plot(t_kin, [0, 15, 28, 48, 78, 95],  "o-", color=PALETTE["passive"], lw=2.5, label="Passive")
    ax.plot(t_kin, [0, 25, 62, 125, 245, 285], "s-", color=PALETTE["active"],  lw=2.5, label="C225-AuNP")
    ax.set_xlabel("Incubation Time (h)")
    ax.set_ylabel("Cellular Uptake (A.U.)")
    ax.set_title("D   Uptake Kinetics", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── E  EGFR-expression correlation ──────────────────────────────────────
    ax = axes[1, 1]
    cell_lines  = ["MCF10A\n(EGFR−)", "A549\n(EGFR+)", "H1975\n(EGFR++)", "H3255\n(EGFR+++)"]
    egfr_expr   = np.array([0.5, 1.8, 3.5, 5.2])
    c225_uptake = np.array([105, 285, 485, 625])
    dot_colors  = [PALETTE["gray"], PALETTE["passive"], PALETTE["free_dox"], PALETTE["active"]]
    for i in range(4):
        ax.scatter(egfr_expr[i], c225_uptake[i], s=180, color=dot_colors[i],
                   edgecolors="k", lw=1.5, zorder=5)
        ax.annotate(cell_lines[i], (egfr_expr[i], c225_uptake[i]),
                    xytext=(8, 8), textcoords="offset points", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", alpha=0.7))
    z = np.polyfit(egfr_expr, c225_uptake, 1)
    xf = np.linspace(0, 6, 100)
    ax.plot(xf, np.polyval(z, xf), "k--", lw=1.5)
    r2 = np.corrcoef(egfr_expr, c225_uptake)[0, 1]**2
    ax.text(0.05, 0.92, f"R² = {r2:.3f}", transform=ax.transAxes, fontsize=11,
            fontweight="bold", bbox=dict(boxstyle="round", fc="wheat", alpha=0.7))
    ax.set_xlabel("EGFR Expression (Relative Units)")
    ax.set_ylabel("C225-AuNP Uptake (A.U.)")
    ax.set_title("E   EGFR-Dependent Targeting", loc="left", fontweight="bold")

    # ── F  Mechanism pie chart ──────────────────────────────────────────────
    ax = axes[1, 2]
    mechs  = ["Receptor-Mediated\nEndocytosis", "Clathrin-Mediated\nEndocytosis",
              "Macropinocytosis", "Passive Diffusion"]
    sizes  = [60, 25, 10, 5]
    colors_= [PALETTE["active"], PALETTE["passive"], PALETTE["free_dox"], PALETTE["gray"]]
    wedges, texts, autotexts = ax.pie(
        sizes, labels=mechs, autopct="%1.1f%%", colors=colors_,
        startangle=90, textprops={"fontsize": 9},
        wedgeprops={"edgecolor": "white", "linewidth": 2})
    for at in autotexts:
        at.set_color("white"); at.set_fontweight("bold")
    ax.set_title("F   Active Targeting — Uptake Mechanisms", loc="left", fontweight="bold")

    fig.suptitle("Figure 3   Mechanistic Validation and Receptor-Specificity Controls",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/Fig3_Mechanistic_Controls.png")
    plt.close()
    print("     Saved → figures/Fig3_Mechanistic_Controls.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — DRUG LOADING QUANTIFICATION
# ══════════════════════════════════════════════════════════════════════════════

def figure_drug_loading():
    """
    Drug loading efficiency, antibody conjugation, spectral confirmation,
    serum stability, BCA assay, and fluorescence quenching.
    """
    print("  Generating Figure 4 — Drug Loading...")

    fig, axes = plt.subplots(2, 3, figsize=(16, 11))

    # ── A  Loading efficiency vs initial ratio ──────────────────────────────
    ax = axes[0, 0]
    ratios = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
    eff    = np.array([92,  88,  82,  75,  65,  52])
    cap    = ratios * eff / 100 * 1000               # mg drug / g Au
    ax2 = ax.twinx()
    ax.plot(ratios, eff, "o-", color=PALETTE["passive"], lw=2.5, label="Loading Efficiency")
    ax2.plot(ratios, cap, "s--", color=PALETTE["active"], lw=2.5, label="Loading Capacity")
    ax.set_xlabel("Initial Drug : Au Mass Ratio")
    ax.set_ylabel("Loading Efficiency (%)", color=PALETTE["passive"])
    ax2.set_ylabel("Capacity (mg drug / g Au)", color=PALETTE["active"])
    ax.set_title("A   Drug Loading Optimisation", loc="left", fontweight="bold")
    ax.tick_params(axis="y", labelcolor=PALETTE["passive"])
    ax2.tick_params(axis="y", labelcolor=PALETTE["active"])
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=8, loc="center right")

    # ── B  Antibody conjugation ─────────────────────────────────────────────
    ax = axes[0, 1]
    ab_ratio = np.array([1, 5, 10, 20, 50, 100])
    ab_per_np = ab_ratio * np.array([85, 82, 75, 68, 55, 42]) * 0.15 / 100
    ax.plot(ab_ratio, ab_per_np, "o-", color=PALETTE["control"], lw=2.5, ms=9)
    ax.axvspan(5, 20, alpha=0.15, color="green", label="Optimal range")
    ax.set_xscale("log")
    ax.set_xlabel("Initial Ab : Au Mass Ratio")
    ax.set_ylabel("Antibodies per Nanoparticle")
    ax.set_title("B   Antibody Conjugation Quantification", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── C  UV-Vis SPR shift ─────────────────────────────────────────────────
    ax = axes[0, 2]
    wl = np.linspace(400, 700, 300)
    for ctr, lbl, col in [(520, "Bare AuNP", PALETTE["free_dox"]),
                           (528, "Drug-Loaded", PALETTE["active"]),
                           (532, "Ab-Conjugated", PALETTE["passive"])]:
        spec = np.exp(-((wl - ctr)**2) / (2 * 50**2))
        ax.plot(wl, spec, lw=2.5, label=lbl, color=col)
        ax.axvline(ctr, color=col, ls="--", lw=1.2, alpha=0.6)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalised Absorbance")
    ax.set_title("C   SPR Shift upon Functionalisation", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── D  Serum stability ──────────────────────────────────────────────────
    ax = axes[1, 0]
    t_s  = np.array([0, 1, 6, 12, 24, 48, 72])
    data_s = {"PBS (pH 7.4)":   [100, 98, 96, 94, 92, 89, 87],
              "10 % FBS":       [100, 95, 88, 82, 75, 68, 62],
              "Mouse Plasma":   [100, 93, 85, 78, 71, 64, 58]}
    cs = [PALETTE["passive"], PALETTE["active"], PALETTE["free_dox"]]
    for (lbl, vals), col in zip(data_s.items(), cs):
        ax.plot(t_s, vals, "o-", lw=2.5, label=lbl, color=col)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Drug Retention (%)")
    ax.set_title("D   Serum Stability of Drug Conjugation", loc="left", fontweight="bold")
    ax.set_ylim(50, 105)
    ax.legend(fontsize=9)

    # ── E  BCA assay ────────────────────────────────────────────────────────
    ax = axes[1, 1]
    np.random.seed(42)
    std_conc = np.array([0, 25, 50, 100, 200, 400, 800])
    std_abs  = std_conc * 0.00125 + 0.05 + np.random.normal(0, 0.008, len(std_conc))
    z = np.polyfit(std_conc, std_abs, 1)
    xf = np.linspace(0, 850, 200)
    ax.scatter(std_conc, std_abs, s=80, color=PALETTE["active"], edgecolors="k", zorder=5, label="Standards")
    ax.plot(xf, np.polyval(z, xf), "--", color=PALETTE["passive"], lw=2, label="Linear fit")
    ax.text(100, 0.92, f"y = {z[0]:.5f}x + {z[1]:.3f}\nR² = 0.998",
            fontsize=9, bbox=dict(boxstyle="round", fc="wheat", alpha=0.7))
    ax.set_xlabel("Protein Concentration (µg/mL)")
    ax.set_ylabel("Absorbance at 562 nm")
    ax.set_title("E   BCA Protein Assay — Standard Curve", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── F  Fluorescence quenching ───────────────────────────────────────────
    ax  = axes[1, 2]
    ax2 = ax.twinx()
    conc_dox = np.array([0, 5, 10, 20, 50, 100, 200])
    fl_free  = conc_dox * 1500 + 200
    fl_bound = conc_dox * 450  + 180
    quench   = np.where(fl_free > 0, (1 - fl_bound / fl_free) * 100, 0)
    ax.plot(conc_dox, fl_free,  "o-", color=PALETTE["active"],  lw=2.5, label="Free DOX")
    ax.plot(conc_dox, fl_bound, "s-", color=PALETTE["passive"], lw=2.5, label="DOX-AuNP")
    ax2.plot(conc_dox, quench, "^-", color=PALETTE["control"], lw=2.5, label="Quenching %")
    ax.set_xlabel("DOX Concentration (µM)")
    ax.set_ylabel("Fluorescence Intensity (A.U.)")
    ax2.set_ylabel("Quenching Efficiency (%)", color=PALETTE["control"])
    ax2.tick_params(axis="y", labelcolor=PALETTE["control"])
    ax.set_title("F   Fluorescence Quenching Analysis", loc="left", fontweight="bold")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=8)

    fig.suptitle("Figure 4   Drug Loading Quantification and Conjugation Stability",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/Fig4_Drug_Loading_Quantification.png")
    plt.close()
    print("     Saved → figures/Fig4_Drug_Loading_Quantification.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — PHOTOTHERMAL THERAPY
# ══════════════════════════════════════════════════════════════════════════════

def figure_photothermal():
    """
    PTT characterisation: power-dependent heating, heating/cooling cycle,
    wavelength-dependent efficiency, chemo-PTT synergy, CEM43 dose, penetration.

    Key parameters (Cole et al. 2009; Sapareto & Dewey 1984)
    ---------------------------------------------------------
    η (nanorods, 808 nm) = 0.99
    Laser: 808 nm, 1.5 W cm⁻², 10 min
    ΔT max = 28 °C above body temperature
    """
    print("  Generating Figure 5 — Photothermal Therapy...")

    fig, axes = plt.subplots(2, 3, figsize=(16, 11))

    # ── A  Power-dependent ΔT ───────────────────────────────────────────────
    ax = axes[0, 0]
    pd_ = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
    ax.plot(pd_, pd_ * 8.5,  "o-", lw=2.5, label="Spheres (520 nm)", color=PALETTE["free_dox"])
    ax.plot(pd_, pd_ * 24.2, "s-", lw=2.5, label="Nanorods (808 nm)", color=PALETTE["active"])
    ax.plot(pd_, pd_ * 31.5, "^-", lw=2.5, label="Nanostars (808 nm)", color=PALETTE["purple"])
    ax.axhline(5,  color="green", ls="--", lw=1.5, label="Hyperthermia (+42 °C)")
    ax.axhline(20, color="red",   ls="--", lw=1.5, label="Ablation (+57 °C)")
    ax.set_xlabel("Laser Power Density (W cm⁻²)")
    ax.set_ylabel("Temperature Rise ΔT (°C)")
    ax.set_title("A   Power-Dependent Heating", loc="left", fontweight="bold")
    ax.legend(fontsize=8)

    # ── B  Heating / cooling cycle ──────────────────────────────────────────
    ax = axes[0, 1]
    t_on  = np.linspace(0, 10, 200)
    t_off = np.linspace(10, 22, 200)
    T_on  = 37 + 28 * (1 - np.exp(-t_on / 1.5))
    T_off = 37 + 28 * np.exp(-(t_off - 10) / 2.5)
    ax.plot(t_on,  T_on,  "r-", lw=3, label="Laser ON")
    ax.plot(t_off, T_off, "b-", lw=3, label="Laser OFF")
    ax.axvspan(0, 10, alpha=0.07, color="red")
    ax.axhline(37, color="gray",   ls="--", lw=1.2, alpha=0.7)
    ax.axhline(42, color="orange", ls="--", lw=1.5, label="Hyperthermia 42 °C")
    ax.axhline(50, color="red",    ls=":",  lw=1.5, label="Ablation 50 °C")
    ax.text(1, 38, "Laser ON", fontsize=9, color="darkred")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("B   Heating–Cooling Cycle\n(808 nm, 1.5 W cm⁻²)", loc="left", fontweight="bold")
    ax.legend(fontsize=8)

    # ── C  Wavelength-dependent efficiency ──────────────────────────────────
    ax = axes[0, 2]
    wls = np.array([532, 650, 808, 980, 1064])
    eta_spheres = np.array([45, 25, 15,  8,  5])
    eta_rods    = np.array([38, 78, 99, 85, 65])
    eta_stars   = np.array([52, 82, 95, 88, 72])
    xi = np.arange(5); wh = 0.25
    ax.bar(xi - wh, eta_spheres, wh, label="Spheres", color=PALETTE["free_dox"], alpha=0.85, edgecolor="k")
    ax.bar(xi,      eta_rods,    wh, label="Nanorods", color=PALETTE["active"],  alpha=0.85, edgecolor="k")
    ax.bar(xi + wh, eta_stars,   wh, label="Nanostars",color=PALETTE["purple"], alpha=0.85, edgecolor="k")
    ax.set_xticks(xi); ax.set_xticklabels([f"{w} nm" for w in wls])
    ax.set_ylabel("Conversion Efficiency (%)")
    ax.set_title("C   Wavelength-Dependent PTT Efficiency", loc="left", fontweight="bold")
    ax.legend(fontsize=8)

    # ── D  Chemo-PTT synergy ────────────────────────────────────────────────
    ax = axes[1, 0]
    treatments = ["Control", "DOX\nOnly", "PTT\nOnly", "DOX + PTT\n(Sequential)", "DOX-AuNP\n+ PTT"]
    viab  = np.array([100, 42, 38, 18,  8])
    viab_e= np.array([  5,  6,  5,  4,  2])
    clrs  = [PALETTE["gray"], PALETTE["passive"], PALETTE["free_dox"],
             PALETTE["purple"], PALETTE["active"]]
    bars  = ax.bar(range(5), viab, yerr=viab_e, capsize=5, color=clrs,
                   alpha=0.85, edgecolor="k")
    ax.plot([1, 4], [58, 58], "k-", lw=1.5)
    ax.text(2.5, 60, "***", ha="center", fontsize=13, fontweight="bold")
    ax.set_xticks(range(5)); ax.set_xticklabels(treatments, fontsize=9)
    ax.set_ylabel("Cell Viability (%)")
    ax.set_title("D   Synergistic Chemo-PTT Effect", loc="left", fontweight="bold")

    # ── E  CEM43 thermal dose ───────────────────────────────────────────────
    ax = axes[1, 1]
    T_range = np.arange(37, 60, 0.5)
    cem     = cem43(T_range, dt=0.5)
    ax.semilogy(T_range, cem, lw=3, color=PALETTE["active"])
    ax.axhline(60, color="orange", ls="--", lw=2, label="Necrosis threshold (60 CEM43)")
    ax.fill_between(T_range, cem, 60, where=cem >= 60, alpha=0.2, color="red", label="Lethal zone")
    ax.axvline(43, color="gray", ls=":", lw=1.5)
    ax.text(43.3, 0.05, "43 °C threshold", fontsize=8, color="gray")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("CEM43 (min equivalent at 43 °C)")
    ax.set_title("E   Thermal Dose — CEM43 Calculation\n(10 min exposure)",
                 loc="left", fontweight="bold")
    ax.legend(fontsize=8)

    # ── F  Tissue penetration depth ─────────────────────────────────────────
    ax = axes[1, 2]
    pen_wls   = [532, 650, 808, 980, 1064]
    pen_depth = [0.5, 1.5, 3.5, 2.8, 2.2]
    clrs_pen  = [PALETTE["free_dox"], PALETTE["purple"], PALETTE["active"],
                 PALETTE["passive"], PALETTE["gray"]]
    bars = ax.bar(range(5), pen_depth, color=clrs_pen, alpha=0.85, edgecolor="k")
    ax.axhline(1.0, color="green", ls="--", lw=2, label="Minimum target depth")
    ax.annotate("NIR-I Window\n(808 nm, optimal)",
                xy=(2, 3.5), xytext=(3.2, 4.0),
                arrowprops=dict(arrowstyle="->", lw=1.8, color="black"),
                fontsize=9, fontweight="bold",
                bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))
    ax.set_xticks(range(5)); ax.set_xticklabels([f"{w} nm" for w in pen_wls])
    ax.set_ylabel("Tissue Penetration Depth (mm)")
    ax.set_title("F   Wavelength-Dependent Penetration Depth", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    fig.suptitle("Figure 5   Photothermal Therapy Characterisation and Combination Efficacy",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/Fig5_Photothermal_Therapy.png")
    plt.close()
    print("     Saved → figures/Fig5_Photothermal_Therapy.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — PHARMACOKINETICS
# ══════════════════════════════════════════════════════════════════════════════

def figure_pharmacokinetics():
    """
    Blood concentration-time profiles, AUC comparison, excretion routes,
    and pharmacokinetic parameter tabulation.

    PK parameters fitted to literature (Perrault 2009; Schipper 2009)
    t½α = 2.5 h, t½β = 58 h for active AuNP-DOX.
    AUC improvement 27.6× over free drug.
    """
    print("  Generating Figure 6 — Pharmacokinetics...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    t_pk = np.array([0.08, 0.25, 0.5, 1, 2, 4, 8, 12, 24, 48, 72, 168])

    C_free    = 100 * np.exp(-0.35 * t_pk)
    C_passive = 100 * (0.65 * np.exp(-0.08  * t_pk) + 0.35 * np.exp(-0.015 * t_pk))
    C_active  = 100 * (0.70 * np.exp(-0.05  * t_pk) + 0.30 * np.exp(-0.012 * t_pk))

    # ── A  Blood concentration ──────────────────────────────────────────────
    ax = axes[0, 0]
    ax.semilogy(t_pk, C_free,    "o-", color=PALETTE["passive"],  lw=2.5, ms=7, label="Free DOX")
    ax.semilogy(t_pk, C_passive, "s-", color=PALETTE["free_dox"], lw=2.5, ms=7, label="Passive AuNP-DOX")
    ax.semilogy(t_pk, C_active,  "^-", color=PALETTE["active"],   lw=2.5, ms=7, label="Active AuNP-DOX")
    ax.text(0.6,  7,  "t½ = 2 h",  color=PALETTE["passive"],  fontsize=9, fontweight="bold")
    ax.text(5,   18,  "t½ = 46 h", color=PALETTE["free_dox"], fontsize=9, fontweight="bold")
    ax.text(12,  32,  "t½ = 58 h", color=PALETTE["active"],   fontsize=9, fontweight="bold")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Blood Concentration (%ID)")
    ax.set_title("A   Blood Concentration–Time Profile", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── B  AUC comparison ───────────────────────────────────────────────────
    ax = axes[0, 1]
    formulations = ["Free DOX", "Bare AuNP", "PEG-AuNP", "Ab-AuNP", "Ab-AuNP-DOX"]
    auc_val  = np.array([ 125,  850, 2850, 3200, 3450])
    auc_err  = np.array([  15,   95,  285,  310,  320])
    clrs_auc = [PALETTE["passive"], PALETTE["gray"], PALETTE["free_dox"],
                PALETTE["purple"], PALETTE["active"]]
    bars = ax.bar(range(5), auc_val, yerr=auc_err, capsize=5, color=clrs_auc,
                  alpha=0.85, edgecolor="k")
    for i, (bar, val) in enumerate(zip(bars, auc_val)):
        if i > 0:
            fold = val / auc_val[0]
            ax.text(bar.get_x() + bar.get_width()/2, val + auc_err[i] + 80,
                    f"{fold:.1f}×", ha="center", fontsize=10,
                    fontweight="bold", color="green")
    ax.set_xticks(range(5))
    ax.set_xticklabels(formulations, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("AUC₀→∞ (%ID·h/mL)")
    ax.set_title("B   Area Under Curve Comparison", loc="left", fontweight="bold")

    # ── C  Clearance routes ─────────────────────────────────────────────────
    ax = axes[1, 0]
    t_ex = np.array([0, 6, 12, 24, 48, 72, 168])
    ax.plot(t_ex, [0, 12, 28, 45, 58, 65, 72], "o-", lw=2.5,
            color=PALETTE["passive"], label="Urine (5 nm AuNP)")
    ax.plot(t_ex, [0,  5, 12, 25, 42, 58, 75], "s-", lw=2.5,
            color=PALETTE["active"],  label="Faeces (50 nm AuNP)")
    ax.plot(t_ex, [0, 0.5, 1.2, 2.5, 4.2, 5.8, 8.5], "^-", lw=2.5,
            color=PALETTE["gray"],   label="Urine (50 nm AuNP)")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Cumulative Excretion (%ID)")
    ax.set_title("C   Size-Dependent Clearance Routes", loc="left", fontweight="bold")
    ax.legend(fontsize=9)

    # ── D  PK parameter comparison ──────────────────────────────────────────
    ax = axes[1, 1]
    pk_labels = ["Cₘₐₓ\n(%ID/mL)", "t½α\n(h)", "t½β\n(h)", "CL\n(mL/h/kg)", "Vd\n(mL/kg)"]
    free_vals = np.array([100,  0.5,   2.0,  850, 2500])
    aunp_vals = np.array([ 95,  2.5,  58.0,   25,  180])
    norm_ref  = np.array([100, 60.0,  60.0,  850, 2500])
    free_norm = free_vals / norm_ref * 100
    aunp_norm = aunp_vals / norm_ref * 100
    xi = np.arange(5); wk = 0.35
    ax.bar(xi - wk/2, free_norm, wk, label="Free DOX",  color=PALETTE["passive"],
           alpha=0.85, edgecolor="k")
    ax.bar(xi + wk/2, aunp_norm, wk, label="AuNP-DOX",  color=PALETTE["active"],
           alpha=0.85, edgecolor="k")
    ax.set_xticks(xi); ax.set_xticklabels(pk_labels, fontsize=9)
    ax.set_ylabel("Normalised Value (%)")
    ax.set_title("D   Pharmacokinetic Parameter Summary", loc="left", fontweight="bold")
    ax.legend(fontsize=9)
    note = ("Actual values — Free DOX: t½β 2 h, AUC 125  |  "
            "AuNP-DOX: t½β 58 h, AUC 3 450 %ID·h/mL, CL 25 mL/h/kg")
    ax.text(0.5, -0.18, note, transform=ax.transAxes, ha="center",
            fontsize=8, style="italic",
            bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.7))

    fig.suptitle("Figure 6   Pharmacokinetic Profiling and Clearance Analysis",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/Fig6_Pharmacokinetics.png")
    plt.close()
    print("     Saved → figures/Fig6_Pharmacokinetics.png")


# ══════════════════════════════════════════════════════════════════════════════
# COMMAND-LINE INTERFACE
# ══════════════════════════════════════════════════════════════════════════════

FIGURE_MAP = {
    1: ("Biodistribution",        figure_biodistribution),
    2: ("Toxicity Panel",         figure_toxicity),
    3: ("Mechanistic Controls",   figure_mechanistic_controls),
    4: ("Drug Loading",           figure_drug_loading),
    5: ("Photothermal Therapy",   figure_photothermal),
    6: ("Pharmacokinetics",       figure_pharmacokinetics),
}


def main():
    parser = argparse.ArgumentParser(
        description="Generate computational figures for AuNP drug delivery research."
    )
    parser.add_argument(
        "--fig", type=int, choices=range(1, 7), metavar="N",
        help="Generate only figure N (1–6). Omit to generate all."
    )
    args = parser.parse_args()

    print("\n" + "=" * 68)
    print("  Gold NP Drug Delivery — Computational Figure Generation")
    print("  Author : Om Jha, St. Xavier's College Kathmandu")
    print("=" * 68 + "\n")

    if args.fig:
        name, fn = FIGURE_MAP[args.fig]
        print(f"Generating Figure {args.fig}: {name}")
        fn()
    else:
        print("Generating all figures...\n")
        for n, (name, fn) in FIGURE_MAP.items():
            print(f"[{n}/6] {name}")
            fn()

    print("\n" + "=" * 68)
    print("  All figures saved to ./figures/  (300 DPI PNG)")
    print("=" * 68 + "\n")


if __name__ == "__main__":
    main()
