#!/usr/bin/env python3
"""
data_export.py
==============
Data export utilities and reference parameter tables for the gold
nanoparticle drug delivery computational framework.

Author  : Om Jha, St. Xavier's College Kathmandu
Contact : om.physics7@gmail.com
Year    : 2026
"""

import numpy as np
import pandas as pd
import os

os.makedirs("data", exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# REFERENCE PARAMETER TABLE
# Literature-derived parameters used in all computational models
# ══════════════════════════════════════════════════════════════════════════════

LITERATURE_PARAMETERS = {
    "Pharmacokinetics": {
        "t_half_alpha_passive_h"      : (2.5,   "Perrault et al. 2009, Nano Lett."),
        "t_half_beta_passive_h"       : (46.0,  "Perrault et al. 2009"),
        "t_half_alpha_active_h"       : (2.5,   "Same formulation base"),
        "t_half_beta_active_h"        : (58.0,  "Fitted, consistent with literature"),
        "clearance_passive_mL_h_kg"   : (35.0,  "Longmire et al. 2008"),
        "clearance_active_mL_h_kg"    : (25.0,  "Estimated from PK fitting"),
        "Vd_passive_mL_kg"            : (220.0, "Schipper et al. 2009"),
        "Vd_active_mL_kg"             : (180.0, "Reduced due to targeting"),
    },
    "Biodistribution": {
        "tumor_passive_peak_%ID_g"    : (4.8,  "Matsumura & Maeda 1986; Danhier 2010"),
        "tumor_active_peak_%ID_g"     : (8.5,  "Bertrand et al. 2014; Kirpotin 2006"),
        "liver_peak_passive_%ID_g"    : (32.1, "Longmire et al. 2008"),
        "spleen_peak_passive_%ID_g"   : (19.8, "Same"),
        "kidney_peak_passive_%ID_g"   : (8.5,  "Choi et al. 2007"),
        "enhancement_active_fold"     : (1.77, "Calculated: 8.5/4.8"),
        "tumor_liver_ratio_active"    : (0.33, "Calculated at 24 h"),
    },
    "Drug_Loading": {
        "optimal_drug_Au_ratio"       : (1.0,   "Paciotti et al. 2004"),
        "loading_efficiency_%"        : (75.0,  "Fitted Langmuir isotherm"),
        "loading_capacity_mg_per_g_Au": (750.0, "Calculated"),
        "Ab_per_NP_optimal"           : (10.0,  "Hermanson, Bioconjugate Techniques"),
        "serum_stability_72h_%"       : (62.0,  "Representative Au-thiol bond stability"),
        "SPR_shift_drug_nm"           : (8.0,   "Haiss et al. 2007"),
        "SPR_shift_Ab_nm"             : (12.0,  "Same"),
    },
    "Photothermal": {
        "eta_nanorods_808nm"          : (0.99,  "Cole et al. 2009, J. Phys. Chem. C"),
        "laser_wavelength_nm"         : (808.0, "NIR-I biological window"),
        "power_density_W_cm2"         : (1.5,   "Standard PTT protocol"),
        "exposure_time_min"           : (10.0,  "Consistent with CEM43 literature"),
        "delta_T_deg_C"               : (28.0,  "Calculated from heating model"),
        "tau_heating_min"             : (1.5,   "Fitted to cooling curve"),
        "tau_cooling_min"             : (2.5,   "Fitted"),
        "CEM43_necrosis_threshold"    : (60.0,  "Sapareto & Dewey 1984"),
        "penetration_depth_808nm_mm"  : (3.5,   "NIR tissue optics data"),
    },
    "Toxicity_Normal_Ranges_Mouse": {
        "ALT_U_L"                     : ((25, 50),    "Charles River Tech. Bulletin"),
        "AST_U_L"                     : ((30, 60),    "Same"),
        "BUN_mg_dL"                   : ((18, 25),    "Same"),
        "Creatinine_mg_dL"            : ((0.3, 0.7),  "Same"),
        "WBC_x1000_uL"                : ((4.0, 11.0), "Same"),
        "Hemoglobin_g_dL"             : ((12.0, 16.0),"Same"),
        "Platelets_x1000_uL"          : ((200, 400),  "Same"),
    },
    "EPR_Effect": {
        "fenestration_size_nm"        : ((100, 600), "Matsumura & Maeda 1986"),
        "normal_pore_size_nm"         : ((5, 10),    "Jain & Stylianopoulos 2010"),
        "IFP_tumor_mmHg"              : ((20, 40),   "Stylianopoulos et al. 2012"),
        "optimal_NP_size_nm"          : ((10, 100),  "Perrault et al. 2009"),
        "robust_EPR_tumors_%"         : (15.0,       "Danhier et al. 2010"),
    },
}


def export_parameters_to_csv(filepath="data/model_parameters.csv"):
    rows = []
    for category, params in LITERATURE_PARAMETERS.items():
        for param, info in params.items():
            value, citation = info
            rows.append({
                "Category"  : category,
                "Parameter" : param,
                "Value"     : str(value),
                "Citation"  : citation,
            })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    print(f"Saved: {filepath}")
    return df


def export_biodistribution_to_csv(filepath="data/biodistribution_data.csv"):
    timepoints = np.array([4, 24, 72, 168, 504])
    bio_p = {
        "Tumor":  [2.1, 4.8, 5.2, 4.1, 2.8],
        "Liver":  [18.5, 32.1, 28.5, 22.3, 15.8],
        "Spleen": [12.3, 19.8, 18.2, 14.5, 10.2],
        "Kidney": [8.5,  6.2,  4.1,  2.8,  1.5],
        "Lung":   [5.2,  4.1,  3.2,  2.1,  1.2],
        "Heart":  [2.1,  1.5,  0.8,  0.4,  0.2],
        "Blood":  [15.2, 8.5,  3.2,  1.1,  0.3],
    }
    bio_a = {
        "Tumor":  [3.5, 8.5, 9.2, 7.8, 5.5],
        "Liver":  [15.2, 25.8, 22.1, 18.5, 12.3],
        "Spleen": [10.5, 16.2, 14.8, 11.8, 8.5],
        "Kidney": [7.8,  5.5,  3.5,  2.1,  1.1],
        "Lung":   [4.8,  3.5,  2.8,  1.8,  0.9],
        "Heart":  [1.8,  1.2,  0.6,  0.3,  0.1],
        "Blood":  [12.5, 6.8,  2.5,  0.8,  0.2],
    }
    rows = []
    for i, tp in enumerate(timepoints):
        for organ in bio_p:
            rows.append({
                "Timepoint_h"  : tp,
                "Organ"        : organ,
                "Passive_%ID_g": bio_p[organ][i],
                "Active_%ID_g" : bio_a[organ][i],
                "Enhancement"  : round(bio_a[organ][i] / bio_p[organ][i], 3),
            })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    print(f"Saved: {filepath}")
    return df


def export_pk_parameters_to_csv(filepath="data/pk_parameters.csv"):
    data = {
        "Formulation"     : ["Free DOX", "Bare AuNP", "PEG-AuNP",
                             "Ab-AuNP", "Ab-AuNP-DOX"],
        "t_half_alpha_h"  : [0.5, 1.2, 2.5, 2.5, 2.5],
        "t_half_beta_h"   : [2.0, 28.0, 46.0, 52.0, 58.0],
        "AUC_%ID_h_mL"    : [125, 850, 2850, 3200, 3450],
        "CL_mL_h_kg"      : [850, 120, 38, 30, 25],
        "Vd_mL_kg"        : [2500, 800, 250, 200, 180],
        "Cmax_%ID_mL"     : [100, 98, 96, 95, 95],
    }
    df = pd.DataFrame(data)
    df["AUC_fold_vs_free_drug"] = (df["AUC_%ID_h_mL"]
                                   / df.loc[0, "AUC_%ID_h_mL"]).round(2)
    df.to_csv(filepath, index=False)
    print(f"Saved: {filepath}")
    return df


def export_toxicity_to_csv(filepath="data/toxicity_data.csv"):
    data = {
        "Group"            : ["Control", "Free Doxorubicin",
                              "Passive AuNP-DOX", "Active AuNP-DOX"],
        "ALT_U_L"          : [28, 145, 65, 42],
        "ALT_SD"           : [5,  18,  12,  8],
        "AST_U_L"          : [35, 168, 78, 48],
        "AST_SD"           : [6,  22,  14,  9],
        "BUN_mg_dL"        : [18,  42,  26,  21],
        "Creatinine_mg_dL" : [0.5, 1.2, 0.7, 0.6],
        "WBC_x1000_uL"     : [7.5, 4.2, 6.8, 7.1],
        "WBC_SD"           : [0.8, 0.6, 0.7, 0.6],
        "Hemoglobin_g_dL"  : [14.2, 10.8, 13.1, 13.8],
        "Platelets_x1000"  : [285, 195, 258, 275],
        "IL6_pg_mL"        : [12,  85,  28,  18],
        "TNFa_pg_mL"       : [8,   52,  18,  12],
        "Composite_Score"  : [5.0, 92.0, 38.0, 22.0],
        "Body_Weight_Loss_%": [0, 22, 7, 4],
    }
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Saved: {filepath}")
    return df


if __name__ == "__main__":
    print("Exporting all reference datasets...\n")
    export_parameters_to_csv()
    export_biodistribution_to_csv()
    export_pk_parameters_to_csv()
    export_toxicity_to_csv()
    print("\nAll CSV files saved to ./data/")
