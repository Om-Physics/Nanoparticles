# Gold Nanoparticles for Targeted Drug Delivery in Lung Cancer

**A Comprehensive Computational Framework for Pharmacokinetic Modelling, Biodistribution Analysis, and Mechanistic Simulation**

**Author:** Om Jha, St. Xavier's College, Kathmandu, Nepal
**Contact:** om.physics7@gmail.com
**ORCID:** 0009-0006-4040-9902
**Version:** 1.0.0
**Year:** 2025

---

## Overview

This repository presents a complete computational framework for analysing gold nanoparticle mediated drug delivery systems in lung adenocarcinoma. The work integrates mathematical pharmacokinetics, biophysical modelling, statistical analysis, and interactive visualisation to provide a rigorous quantitative evaluation of both passive and active tumour targeting strategies.

Lung cancer remains the leading cause of cancer mortality worldwide, accounting for approximately 2.2 million new diagnoses and 1.8 million deaths annually. Conventional chemotherapy is severely limited by poor tumour selectivity and systemic toxicity that significantly diminishes quality of life. Gold nanoparticles offer a compelling platform to overcome these limitations through their exceptional biocompatibility, size-tunable surface plasmon resonance, straightforward surface chemistry, and capacity to integrate therapeutic, diagnostic, and photothermal functions within a single nanoscale construct.

This computational study systematically evaluates the key determinants of nanoparticle performance including organ biodistribution, pharmacokinetic profiles, cellular uptake mechanisms, pH-responsive drug release, photothermal therapy characterisation, and comprehensive toxicity assessment. All quantitative results represent computational projections derived from established models with parameters sourced from peer-reviewed literature. They are intended to support experimental design and provide mechanistic insight rather than to substitute for empirical measurement.

---

## Scientific Background

### Gold Nanoparticle Properties

Gold nanoparticles exhibit surface plasmon resonance arising from the collective oscillation of conduction band electrons driven by incident electromagnetic radiation. For spherical particles approximately 50 nanometres in diameter, this resonance occurs near 520 nanometres in aqueous solution, producing an intense absorption and scattering signature directly observable by ultraviolet-visible spectroscopy. The resonance wavelength red-shifts with increasing particle size and is highly sensitive to surface modifications, providing a convenient spectroscopic handle for monitoring conjugation chemistry. Gold nanorods display two resonance peaks corresponding to transverse and longitudinal plasmon modes, with the longitudinal mode tuneable into the near-infrared biological window between 700 and 1100 nanometres, making them especially attractive for photothermal therapy.

The high surface area to volume ratio of gold nanoparticles permits dense loading of therapeutic cargo, targeting ligands, and stabilising polymers. PEGylation through thiol-terminated polyethylene glycol chains provides a hydrophilic corona that suppresses opsonisation, reduces reticuloendothelial system recognition, and extends blood circulation half-life from minutes for bare particles to 46 hours or more. Cetuximab, a chimeric monoclonal antibody against epidermal growth factor receptor, can be conjugated to the PEG-AuNP surface through established carbodiimide or maleimide coupling chemistry, introducing active targeting capability while preserving colloidal stability.

### The Enhanced Permeability and Retention Effect

The Enhanced Permeability and Retention effect, first described by Matsumura and Maeda in 1986, is the cornerstone of passive tumour targeting. Rapidly proliferating tumour tissue generates a structurally defective vasculature with interendothelial fenestrations ranging from 100 to 600 nanometres, compared to 5 to 10 nanometres in normal capillaries. Nanoparticles within the optimal size window of 10 to 100 nanometres selectively extravasate through these large openings and accumulate in the tumour interstitium. Impaired lymphatic drainage further prevents clearance, sustaining elevated intratumoural concentrations over days. Computational modelling of this process employs first-order differential equations balancing EPR-mediated influx driven by plasma concentration against diffusion-limited efflux.

### Active Targeting via EGFR

Epidermal growth factor receptor is overexpressed in 40 to 80 percent of non-small cell lung cancers, making it a clinically validated target for antibody-conjugated nanoparticles. Cetuximab binding to EGFR domain III triggers receptor-mediated clathrin-dependent endocytosis, internalising the nanoparticle into early endosomes that progressively acidify to pH 5.0 through proton pump activity. Doxorubicin conjugated through pH-sensitive hydrazone bonds is selectively released in this acidic environment, bypassing P-glycoprotein efflux pumps and delivering drug directly to the nucleus. This active targeting mechanism increases tumour accumulation from approximately 4.8 to 8.5 percent injected dose per gram compared to passive EPR-only delivery, a 77 percent enhancement that translates to substantially reduced systemic exposure and improved therapeutic index.

---

## Repository Structure

```
Nanoparticles/
├── README.md                        This document
├── LICENSE                          MIT License
├── CHANGELOG.md                     Version history and release notes
├── requirements.txt                 Python package dependencies
├── setup.py                         Package installation configuration
├── .gitignore                       Version control exclusions
│
├── src/                             Source code modules
│   ├── __init__.py
│   ├── models.py                    All mathematical and pharmacokinetic models
│   ├── generate_figures.py          Publication figure generation (six figures)
│   ├── statistics_utils.py          Statistical analysis and reporting utilities
│   └── data_export.py               CSV export and literature parameter reference
│
├── docs/                            Scientific documentation
│   ├── RESEARCH_PAPER.md            Complete manuscript with abstract through references
│   ├── CODE_GUIDE.md                Function-level technical documentation
│   ├── FIGURES_GUIDE.md             Panel descriptions, results, and reproduction guide
│   ├── SIMULATION_GUIDE.md          Interactive simulation architecture and usage
│   └── INDEX.md                     Master navigation document
│
├── simulation/
│   └── Gold_NP_Simulation.html      Self-contained interactive mechanism visualisation
│
├── examples/
│   └── quick_start.py               Minimal working example of all core analyses
│
├── tests/
│   └── test_models.py               Unit tests for all mathematical model functions
│
├── figures/                         Generated output directory (300 DPI PNG)
└── data/                            Generated CSV reference datasets
```

---

## Installation

### System Requirements

Python 3.8 or higher is required. The framework has been tested on Windows 10, macOS 11 Big Sur, and Ubuntu 20.04 LTS. A minimum of 4 GB RAM is recommended for figure generation; all analyses complete within seconds on a standard laptop.

### Setting Up the Environment

It is strongly recommended to create an isolated virtual environment before installation to prevent dependency conflicts with other Python projects.

```bash
git clone https://github.com/Om-Physics/Nanoparticles.git
cd Nanoparticles
python -m venv venv
source venv/bin/activate          # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

For development use, install the package in editable mode so that changes to source files take effect immediately without reinstallation.

```bash
pip install -e .
```

Verify correct installation by running the unit test suite.

```bash
pytest tests/ -v
```

A passing result confirms that all mathematical models return biologically plausible values and that the statistical utilities are functioning correctly.

---

## Quick Start

### Generating All Six Publication Figures

```bash
python src/generate_figures.py
```

This command generates all figures simultaneously and saves them as 300 DPI PNG files in the `figures/` directory. To generate a specific figure by number, use the `--fig` flag.

```bash
python src/generate_figures.py --fig 1   # Biodistribution only
python src/generate_figures.py --fig 5   # Photothermal therapy only
```

### Running the Interactive Simulation

Open `simulation/Gold_NP_Simulation.html` directly in any modern web browser. No server, internet connection, or additional software is required. The simulation loads immediately and renders eighteen gold nanoparticles and three A549 lung adenocarcinoma cells through all eight mechanistic phases of the drug delivery process.

### Exporting Reference Datasets

```bash
python src/data_export.py
```

This exports four CSV files to the `data/` directory covering literature model parameters, complete biodistribution datasets, pharmacokinetic parameter comparisons, and toxicity biomarker reference values.

### Running the Worked Example

```bash
python examples/quick_start.py
```

This script demonstrates pharmacokinetic profiling, EPR tumour accumulation, pH-responsive drug release, photothermal CEM43 dose calculation, IC50 determination, and comprehensive statistical analysis, printing formatted results to the terminal.

---

## Computational Modules

### Mathematical Models (`src/models.py`)

The models module implements all governing equations used across the computational analyses. The two-compartment pharmacokinetic model describes blood concentration as a biexponential function with a rapid distribution phase and a slower elimination phase, fitted to PEGylated gold nanoparticle data from Perrault et al. (2009). Tumour accumulation follows a first-order differential equation driven by EPR-mediated influx from the plasma compartment, numerically integrated using `scipy.odeint`.

Drug release from nanoparticle surfaces is modelled using three complementary approaches: the Korsmeyer-Peppas power law for characterising diffusion mechanism through the release exponent, pH-dependent first-order kinetics with a Hill-type sensitivity function for endosomal release selectivity, and Langmuir adsorption isotherms for surface loading equilibria. Photothermal heating follows an exponential rise model parameterised by conversion efficiency, laser power density, and thermal time constant, with the CEM43 cumulative thermal dose calculated using the Sapareto and Dewey formulation.

Optical and structural characterisation models include Mie quasi-static extinction theory for nanoparticle concentration determination, the Stokes-Einstein equation for hydrodynamic radius from diffusion measurements, and the Scherrer equation for XRD crystallite size analysis. The Starling fluid flux equation quantifies hydrostatic and oncotic pressure contributions to vascular permeability governing nanoparticle extravasation.

### Figure Generation (`src/generate_figures.py`)

The figure generation module serves as the primary entry point for creating all six publication figures. Each figure function is self-contained and independently callable. A consistent global style is applied through `matplotlib.rcParams` using serif fonts, clean axis spines, and a controlled colour palette that distinguishes passive targeting (blue), active targeting (red), free drug (orange), control (green), and supplementary groups (purple, grey).

Figure 1 presents quantitative organ biodistribution across tumour, liver, spleen, kidney, lung, heart, and blood compartments at five time points spanning 4 hours to 21 days, comparing passive PEG-AuNP-DOX and active cetuximab-AuNP-DOX formulations. Figure 2 provides an eight-panel toxicity assessment covering liver enzymes, kidney markers, haematology, inflammatory cytokines, histopathological scoring, body weight, and a composite safety score. Figure 3 validates active targeting specificity through isotype controls, receptor blocking competition, endocytosis pathway inhibitors, time-dependent uptake kinetics, EGFR expression correlation, and pathway mechanism distribution. Figure 4 characterises drug loading through Langmuir efficiency curves, antibody conjugation quantification, SPR spectral shifts, serum stability profiling, BCA calibration, and fluorescence quenching analysis. Figure 5 covers photothermal therapy through laser parameter optimisation, heating-cooling cycle dynamics, wavelength-dependent tissue penetration, chemo-PTT synergy, CEM43 thermal dose mapping, and wavelength comparison across the optical window. Figure 6 presents pharmacokinetic blood concentration-time profiles, AUC fold-improvement comparisons, size-dependent excretion route analysis, and a normalised parameter summary table.

### Statistical Utilities (`src/statistics_utils.py`)

The statistics module provides all significance testing and reporting functions consistent with biomedical publication standards. The two-sample t-test function returns a complete reporting dictionary including the test statistic, degrees of freedom, exact p-value, Cohen's d effect size with verbal interpretation, 95% confidence interval for the difference in means, and a significance label in the standard asterisk notation. One-way ANOVA is supplemented with eta-squared effect size and is designed to accept any number of groups through variable positional arguments.

Additional functions include power analysis for minimum sample size determination given a specified effect size, type I error rate, and target power; Shapiro-Wilk normality testing with interpretation; the Mantel-Haenszel log-rank test for survival curve comparison; Bonferroni correction for multiple comparison families; and comprehensive descriptive statistics reporting mean, median, standard deviation, standard error, and 95% confidence interval.

### Data Export (`src/data_export.py`)

The data export module provides programmatic access to all reference parameter values compiled from the primary literature. The `LITERATURE_PARAMETERS` dictionary organises parameters by category covering pharmacokinetics, biodistribution, drug loading, photothermal characterisation, normal mouse laboratory reference ranges, and EPR effect quantitative descriptors. Each entry pairs the numerical value with the full bibliographic citation enabling traceability back to source publications.

Four export functions generate structured CSV files for all major datasets: the complete parameter reference table, organ biodistribution data at all five time points for both formulations, pharmacokinetic parameter comparisons across five formulations, and toxicity biomarker profiles with standard deviations for all treatment groups.

---

## Key Results Summary

The computational analyses establish several quantitatively significant findings that motivate experimental investigation of the proposed system. Active cetuximab-conjugated AuNP-DOX achieves a peak tumour accumulation of 8.5 percent injected dose per gram at 24 hours compared to 4.8 percent for passive EPR delivery, representing a 77 percent enhancement in tumour targeting efficiency. This comes alongside a meaningful reduction in liver and spleen sequestration of 20 to 25 percent relative to passive formulations, improving the tumour-to-reticuloendothelial organ ratio.

Pharmacokinetic profiling demonstrates a terminal half-life of 58 hours for active AuNP-DOX compared to 2 hours for free doxorubicin, producing an AUC improvement of 27.6-fold. This prolonged circulation directly enables EPR accumulation and active receptor engagement. The composite toxicity score for active AuNP-DOX of 22 out of 100 compares favourably with 92 for free doxorubicin, driven primarily by the substantially reduced liver enzyme elevations and haematological suppression.

pH-responsive doxorubicin release is selective by at least 3-fold between endosomal pH 5.0 and physiological pH 7.4, with 24-hour release of approximately 78 percent at endosomal pH versus less than 18 percent at blood pH. Combination chemo-photothermal therapy using 808 nm NIR irradiation at 1.5 W/cm² for 10 minutes achieves a temperature rise of 28 degrees Celsius and a CEM43 thermal dose exceeding the 60 minute equivalent necrosis threshold, reducing cell viability to 8 percent compared to 38 percent for photothermal treatment alone.

---

## Mathematical Model Reference

| Model | Governing Equation | Source Reference |
|---|---|---|
| Two-compartment pharmacokinetics | C(t) = A·exp(−αt) + B·exp(−βt) | Gibaldi and Perrier, 1982 |
| EPR tumour accumulation | dC/dt = k_in·C_plasma − k_out·C | Cabral et al., 2011 |
| Korsmeyer-Peppas drug release | M_t/M_inf = k·t^n | Korsmeyer et al., 1983 |
| pH-dependent release | k(pH) = k_max / (1 + (pH/pH_ref)^hill) | Framework-derived |
| Four-parameter logistic | Response = bottom + (top−bottom) / (1+(dose/IC50)^hill) | DeLean et al., 1978 |
| Photothermal heating | ΔT(t) = η·I·(1−exp(−t/τ)) | Cole et al., 2009 |
| CEM43 thermal dose | Σ R^(43−T)·Δt | Sapareto and Dewey, 1984 |
| Langmuir adsorption | q = C_max·C / (K_d + C) | Langmuir, 1918 |
| Mie extinction (quasi-static) | ε ∝ R³·ε_m^(3/2)·ε_i / ((ε_r+2ε_m)² + ε_i²) | Haiss et al., 2007 |
| Stokes-Einstein radius | R_h = k_B·T / (6π·η·D) | Einstein, 1905 |
| Scherrer crystallite size | D = K·λ / (β·cosθ) | Scherrer, 1918 |
| Starling fluid flux | J_v = L_p·S·[(P_c−P_i) − σ(π_c−π_i)] | Stylianopoulos et al., 2012 |

---

## Testing

The test suite in `tests/test_models.py` provides 40 unit tests organised across 10 test classes covering pharmacokinetic models, EPR accumulation, drug release kinetics, dose-response fitting, photothermal models, optical and structural models, the EPR efficiency heuristic, Starling fluid flux, statistical utilities, and a mass balance validation. Tests verify mathematical correctness by checking known analytical solutions, biological plausibility by asserting literature-range outputs, and edge case robustness including zero inputs, very large values, and saturation limits.

Run the full suite with verbose output.

```bash
pytest tests/ -v
```

Run a specific test class.

```bash
pytest tests/test_models.py::TestPhotothermal -v
```

---

## Citation

When using this framework in published work, please include the following statement in your methods section and cite the repository.

"Computational pharmacokinetic, biodistribution, toxicity, and mechanistic analyses were performed using a custom Python framework implementing two-compartment pharmacokinetic models, EPR accumulation differential equations, Korsmeyer-Peppas drug release kinetics, and CEM43 photothermal dose calculations. All parameters were derived from peer-reviewed literature on PEGylated and antibody-conjugated gold nanoparticle systems. All results represent computational projections and not experimental measurements. Code and full parameter documentation are available at github.com/Om-Physics/Nanoparticles."

---

## License

This project is distributed under the MIT License. You are free to use, modify, and distribute the code with appropriate attribution. See the LICENSE file for the complete terms. No warranties are provided regarding the accuracy, completeness, or fitness for any particular purpose of the computational outputs.

---

## Contact and Acknowledgments

**Om Jha** — om.physics7@gmail.com — ORCID: 0009-0006-4040-9902
St. Xavier's College, Kathmandu, Nepal

Bug reports and feature requests: github.com/Om-Physics/Nanoparticles/issues

This work builds on the foundational contributions of Matsumura and Maeda in describing the Enhanced Permeability and Retention effect, the extensive experimental literature on PEGylated gold nanoparticle pharmacokinetics, and the open-source scientific computing community whose tools make this level of quantitative analysis accessible to researchers worldwide.

---

*Version 1.0.0 — February 2026*
