# Code Guide

## Technical Documentation for the Gold Nanoparticle Drug Delivery Computational Framework

**Repository:** github.com/Om-Physics/Nanoparticles
**Author:** Om Jha, St. Xavier's College, Kathmandu, Nepal
**Version:** 1.0.0 — February 2026

---

## Overview

This guide provides function-level technical documentation for all four source modules in the framework. Each section describes the mathematical basis, input parameters with units, return values, usage examples, and the primary literature from which equations and parameter values are derived. Readers familiar with the scientific background may navigate directly to the module of interest using the section headings below.

The framework is organised as follows. The `models` module implements every governing equation used in the analyses and is imported by all other modules. The `generate_figures` module calls model functions and produces all six publication figures. The `statistics_utils` module provides standalone significance testing and reporting functions. The `data_export` module supplies reference parameter tables and exports structured datasets to CSV.

---

## Module 1: Mathematical Models (`src/models.py`)

### Pharmacokinetic Functions

**`two_compartment_pk(t, A, alpha, B, beta)`**

Implements the classical two-compartment open pharmacokinetic model describing drug concentration in the central (blood) compartment following intravenous bolus administration.

The governing equation is C(t) = A·exp(−αt) + B·exp(−βt), where the first term represents the rapid distribution phase as drug partitions into peripheral tissues and the second term represents the slower elimination phase governed by renal and hepatic clearance. At time zero, C(0) equals A plus B, which corresponds to the initial concentration in the central compartment immediately after injection.

Parameters A and alpha describe the distribution phase: A is the intercept coefficient in percent injected dose per millilitre, and alpha is the distribution rate constant in reciprocal hours. Parameters B and beta describe the elimination phase with equivalent units. The terminal half-life is computed as ln(2) divided by beta.

For PEGylated AuNP-DOX, fitted values from Perrault et al. (2009) give A of 65, alpha of 0.277 per hour, B of 35, and beta of 0.015 per hour, yielding a terminal half-life of 46 hours. For active C225-AuNP-DOX, beta of 0.012 per hour produces a 58-hour terminal half-life.

```python
import numpy as np
from src.models import two_compartment_pk

t = np.linspace(0, 168, 2000)
C_active = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
```

**`epr_tumor_accumulation(t, k_in, k_out, plasma_params)`**

Solves the first-order differential equation governing EPR-mediated tumour accumulation.

dC_tumour/dt = k_in · C_plasma(t) − k_out · C_tumour

The plasma concentration C_plasma is computed internally from the two-compartment model using the parameter tuple provided in `plasma_params` as (A, alpha, B, beta). The ODE is integrated numerically from time zero using SciPy's `odeint` solver with an initial tumour concentration of zero.

The influx rate constant k_in in reciprocal hours incorporates vascular permeability, available fenestration surface area, and the probability of nanoparticle passage through any given fenestration at a given particle size. The efflux rate constant k_out reflects tumour interstitial clearance and re-entry into lymphatics. For passive formulations, k_in of 0.085 and k_out of 0.015 per hour produce a peak accumulation of 4.8 percent ID per gram at approximately 24 hours. For active formulations, k_in of 0.110 reflects additional receptor-mediated retention.

```python
from src.models import epr_tumor_accumulation

plasma = (70, 0.277, 30, 0.012)
C_tumour = epr_tumor_accumulation(t, k_in=0.110, k_out=0.012, plasma_params=plasma)
peak_val = C_tumour.max()
peak_time = t[C_tumour.argmax()]
```

**`calculate_pk_parameters(t, C)`**

Derives the four primary pharmacokinetic metrics from a concentration-time array: Cmax (peak concentration), Tmax (time of peak), AUC (area under the curve from zero to the last measured time point by trapezoidal integration), and t_half (terminal half-life estimated by log-linear regression of the last three data points). Returns a dictionary with keys Cmax, Tmax, AUC, and t_half.

**`pbpk_organ(t, Q, V, P, plasma_fn)`**

A simplified physiologically-based pharmacokinetic model for organ-specific distribution.

dC_organ/dt = (Q/V) · (C_blood − C_organ/P)

Organ blood flow Q in millilitres per hour, organ volume V in millilitres, and the tissue-plasma partition coefficient P are organ-specific parameters. The function accepts a callable `plasma_fn` that returns the plasma concentration at any time point, enabling integration with the two-compartment plasma model. This model is used for secondary organ distributions but is not the primary model for liver and spleen, which are dominated by RES uptake and require empirical scaling factors.

---

### Drug Release Functions

**`first_order_release(t, k, M_total=100.0)`**

First-order release kinetics: M(t) = M_total · (1 − exp(−k·t)).

This is the simplest release model and describes exponential approach to equilibrium. It is most appropriate for drug dissolved homogeneously within a matrix with no diffusional barriers. The rate constant k in reciprocal hours determines the speed of release; the function returns percent cumulative release clipped at M_total.

**`korsmeyer_peppas_release(t, k, n, M_total=100.0)`**

Korsmeyer-Peppas power law: M_t/M_inf = k·t^n.

The release exponent n provides mechanistic insight into the transport mechanism. For spherical particles, n below 0.43 indicates Fickian diffusion dominated by concentration gradient; n between 0.43 and 0.85 indicates anomalous transport combining diffusion and polymer relaxation; n above 0.85 indicates Case-II transport controlled entirely by polymer chain relaxation rather than diffusion. Fitted values for the AuNP-DOX system at pH 7.4 give k of 0.018 and n of 0.52, indicating anomalous transport.

**`ph_dependent_release(t, pH, k_max=0.05, pH_ref=6.0, hill=3.0, M_total=100.0)`**

Combines a Hill-type pH sensitivity function with first-order release kinetics to model hydrazone bond cleavage.

k(pH) = k_max / (1 + (pH/pH_ref)^hill)

At pH 7.4, the denominator is approximately 5.1 with the default hill coefficient of 3, yielding k(7.4) of approximately 0.010 per hour. At pH 5.0, the denominator approaches 1, yielding k(5.0) approaching k_max of 0.05 per hour. This 5-fold difference in effective rate constant produces the greater than 3-fold selectivity in 24-hour cumulative release. The hill coefficient of 3 provides a sharp sigmoidal transition consistent with the cooperative pH response of hydrazone-based systems.

```python
import numpy as np
from src.models import ph_dependent_release

t = np.linspace(0, 72, 720)
release_blood   = ph_dependent_release(t, pH=7.4)
release_tumour  = ph_dependent_release(t, pH=6.5)
release_endo    = ph_dependent_release(t, pH=5.0)
selectivity = release_endo[-1] / release_blood[-1]
print(f"Endosomal to blood selectivity at 72 h: {selectivity:.1f}x")
```

---

### Dose-Response Functions

**`four_parameter_logistic(dose, bottom, top, IC50, hill)`**

The four-parameter logistic model is the standard for IC50 determination in cytotoxicity assays.

Response = bottom + (top − bottom) / (1 + (dose/IC50)^hill)

The bottom parameter is the minimum response (viability at saturating drug concentration, typically approaching 0 percent). The top parameter is the maximum response at zero drug concentration, typically 100 percent. IC50 is the dose producing half-maximal effect. The hill coefficient, negative for inhibitory responses, determines the steepness of the sigmoidal transition.

At dose equal to IC50, the formula reduces to bottom + (top − bottom)/2, confirming that the response at IC50 is always precisely the midpoint between bottom and top regardless of the hill coefficient.

**`fit_ic50(doses, viabilities)`**

Wraps `four_parameter_logistic` in SciPy's `curve_fit` nonlinear least-squares optimiser and returns a four-element tuple: (IC50, hill, popt_full, r_squared). Initial parameter guesses are derived automatically from the data. Returns R-squared as a goodness-of-fit measure computed from residual and total sums of squares.

---

### Photothermal Functions

**`photothermal_heating(t, eta, I, tau, T_body=37.0)`**

Models temperature rise during continuous-wave laser irradiation.

T(t) = T_body + η · I · (1 − exp(−t/τ))

The photothermal conversion efficiency eta for gold nanorods at 808 nm is 0.99 based on Cole et al. (2009), meaning that 99 percent of absorbed photon energy is converted to heat rather than re-emitted as fluorescence. The power density I in watts per square centimetre and the thermal time constant tau in seconds jointly determine the heating rate and steady-state temperature. At standard conditions (I of 1.5 W/cm², tau of 90 s), steady-state ΔT is 28 degrees Celsius above the baseline body temperature of 37 degrees Celsius.

**`photothermal_cooling(t, T_peak, tau_cool, T_body=37.0)`**

Newton's law of cooling following laser cessation.

T(t) = T_body + (T_peak − T_body) · exp(−t/τ_cool)

The cooling time constant tau_cool of 150 seconds for the standard protocol returns temperature within 0.5 degrees of baseline within 15 minutes of laser shutdown, confirming spatially confined thermal damage.

**`cem43_dose(T_profile, dt=1.0)`**

Sapareto-Dewey cumulative equivalent minutes at 43 degrees Celsius.

CEM43 = Σ R^(43−T) · Δt

where R equals 0.25 below 43 degrees Celsius and 0.5 at or above 43 degrees Celsius. This formulation captures the temperature-dependent rate of protein denaturation and cell killing, normalising diverse time-temperature exposure histories to a single equivalent quantity. The threshold of CEM43 greater than 60 is established as the minimal lethal dose for irreversible tissue damage across multiple tissue types.

```python
import numpy as np
from src.models import photothermal_heating, cem43_dose

t_s = np.linspace(0, 600, 600)
T   = photothermal_heating(t_s, eta=0.99, I=1.5, tau=90, T_body=37.0)
cem = cem43_dose(T, dt=1.0)
print(f"Peak temperature: {T.max():.1f} C")
print(f"CEM43 at end of exposure: {cem[-1]:.1f}")
print(f"Necrosis threshold exceeded: {cem[-1] > 60}")
```

---

### Optical and Structural Functions

**`mie_extinction_coefficient(R, epsilon_m, epsilon_r, epsilon_i, wavelength)`**

Quasi-static Mie extinction coefficient for spherical nanoparticles in the dipole approximation, valid when the particle radius is much smaller than the wavelength. Returns the molar extinction coefficient in M⁻¹ cm⁻¹, which for 50 nm gold nanoparticles at 520 nm is typically in the range of 10⁹ to 10¹⁰ M⁻¹ cm⁻¹. Used for nanoparticle concentration determination from UV-Vis absorbance measurements.

**`stokes_einstein_radius(D, T=298.15, eta=0.001)`**

Converts a measured diffusion coefficient in m²/s (from dynamic light scattering or NMR diffusometry) to a hydrodynamic radius in nanometres using the Stokes-Einstein relation. The default solvent viscosity corresponds to water at 25 degrees Celsius. At physiological temperature of 37 degrees Celsius, water viscosity is approximately 0.00069 Pa·s.

**`scherrer_crystallite_size(beta_rad, theta_rad, K=0.9, wavelength_nm=0.154)`**

XRD crystallite size from the Scherrer equation using Cu Kα radiation at 0.154 nm by default. The shape factor K of 0.9 is appropriate for spherical crystallites. Both the peak width beta and Bragg angle theta must be provided in radians. The function returns crystallite size in nanometres.

---

### Adsorption and Transport Functions

**`langmuir_adsorption(C, C_max, K_d)`**

Langmuir saturation isotherm for drug loading on nanoparticle surfaces at equilibrium. At concentrations much less than K_d, loading is linearly proportional to free drug concentration. At concentrations much greater than K_d, loading approaches C_max asymptotically. For the AuNP-DOX system, K_d of 8.2 µg/mL and C_max of 750 µg/mg Au are obtained from isotherm fitting.

**`epr_accumulation_efficiency(diameter_nm, zeta_mV, t_circ_h, K_perm)`**

A composite heuristic function combining four determinants of EPR efficiency: a Gaussian size-selectivity factor centred at 50 nm, a charge stability factor based on zeta potential magnitude, a circulation time factor with a 20-hour characteristic scale, and a permeability coefficient factor. Returns a dimensionless efficiency value between 0 and 1. This function is intended for relative comparison of nanoparticle designs rather than absolute prediction of accumulation values.

**`starling_fluid_flux(L_p, S, P_c, P_i, sigma, pi_c, pi_i)`**

Net fluid flux across capillary walls from the Starling equation, balancing hydrostatic and oncotic pressure driving forces against the osmotic reflection coefficient. In tumours, elevated interstitial fluid pressure (P_i of 20 to 40 mmHg compared to −3 mmHg in normal tissue) reduces the effective hydrostatic pressure gradient and can oppose nanoparticle extravasation despite large fenestration size.

---

## Module 2: Figure Generation (`src/generate_figures.py`)

### Global Configuration

The module sets consistent figure aesthetics through `matplotlib.rcParams` at import time using a serif font family (Times New Roman preferred, DejaVu Serif fallback), axis label size of 12 points, tick label size of 10 points, line width of 2.2 points, and a grid with alpha 0.3 and dashed style. Top and right spines are removed for a clean presentation. All figures are saved at 300 DPI using tight bounding box.

The `PALETTE` dictionary defines six named colours used consistently across all figures: passive targeting in blue (#3498db), active targeting in red (#e74c3c), control in green (#2ecc71), free drug in orange (#f39c12), grey for baseline or unlabelled groups (#95a5a6), and purple for additional treatment arms (#9b59b6).

### Adding a New Figure

To add a new figure to the framework, create a function named `figure_yourtopic()` following the pattern of existing functions. The function should print a progress message, create a figure with descriptive panel labels in bold using the `loc="left"` title alignment convention, apply the global PALETTE, save to `figures/FigN_YourTopic.png` at 300 DPI, close the figure, and print the save path. Add the function to the `FIGURE_MAP` dictionary at the bottom of the module with a sequential integer key.

```python
FIGURE_MAP = {
    ...
    7: ("My New Analysis", figure_mynewanalysis),
}
```

### Command-Line Interface

The `main()` function parses a single optional argument `--fig N` where N is an integer between 1 and 6. When called without arguments, all figures are generated in sequence. The interface is accessible via both `python src/generate_figures.py` and the `generate-figures` console entry point defined in `setup.py`.

---

## Module 3: Statistical Utilities (`src/statistics_utils.py`)

### Function Summary

**`two_sample_ttest(group1, group2, alternative="two-sided", alpha=0.05)`**

Performs an unpaired two-sample t-test (Welch correction is not applied; both groups are assumed to have equal or similar variances, appropriate for matched experimental designs). Returns a dictionary containing: t_stat (four decimal places), p_value (six decimal places), df (integer degrees of freedom), significant (boolean), sig_label (asterisk notation), cohens_d (effect size), effect_size (verbal label: negligible, small, medium, or large), ci_95 (tuple of the 95 percent confidence interval for the difference in means), and mean_diff.

**`one_way_anova(*groups, alpha=0.05)`**

Accepts any number of group arrays as positional arguments. Returns F_stat, p_value, significant, eta_squared, and n_groups. Eta-squared is computed as the ratio of between-group sum of squares to total sum of squares and is reported with four decimal places.

**`sample_size_calculation(effect_size=0.8, alpha=0.05, power=0.80)`**

Iterative minimum sample size search using the non-central t-distribution. Increments n from 2 until achieved power equals or exceeds the specified target. Returns n_per_group, total_n, actual_power, effect_label, and effect_size. Accepts Cohen's d values where 0.2 is small, 0.5 is medium, and 0.8 is large.

**`normality_test(data, alpha=0.05)`**

Shapiro-Wilk test via `scipy.stats.shapiro`. Returns statistic, p_value, is_normal (boolean), and an interpretation string. Most appropriate for sample sizes between 3 and 50; for larger samples the test becomes extremely sensitive to minor deviations and the result should be interpreted with caution.

**`log_rank_test(time1, event1, time2, event2)`**

Mantel-Haenszel log-rank statistic for comparing two survival curves. Event indicators should be 1 for an observed event and 0 for a censored observation. Returns chi2, p_value, significant, and hazard_ratio estimated from event rates.

**`bonferroni_correction(p_values, alpha=0.05)`**

Divides the family-wise error rate alpha by the number of comparisons m to obtain a corrected per-comparison threshold. Also returns adjusted p-values (each raw p-value multiplied by m, clipped at 1) and a boolean significance mask.

**`summary_statistics(data, label="")`**

Comprehensive descriptive statistics: n, mean, median, standard deviation, standard error of the mean, 95% confidence interval using the t-distribution, minimum, and maximum. All numeric values are rounded to four decimal places.

---

## Module 4: Data Export (`src/data_export.py`)

### The LITERATURE_PARAMETERS Dictionary

This dictionary is the authoritative parameter reference for the framework. It is organised into six top-level categories. Every entry is a two-element tuple of (value, citation_string) where value may be a scalar or a range tuple.

The Pharmacokinetics category contains distribution and elimination rate constants, clearance rates, and volumes of distribution for all five formulations. The Biodistribution category contains peak percent ID per gram values for tumour and major organs. The Drug_Loading category contains loading efficiency, conjugation parameters, spectral shift values, and stability metrics. The Photothermal category contains conversion efficiency, laser parameters, heating and cooling time constants, and the CEM43 necrosis threshold. The Toxicity_Normal_Ranges_Mouse category contains reference intervals for ALT, AST, BUN, creatinine, WBC, haemoglobin, and platelets in C57BL/6 mice. The EPR_Effect category contains fenestration size ranges, normal pore sizes, interstitial fluid pressure, optimal size range, and the percentage of tumours exhibiting robust EPR response.

### Export Functions

`export_parameters_to_csv()` flattens the nested parameter dictionary into a tabular CSV with columns Category, Parameter, Value, and Citation.

`export_biodistribution_to_csv()` generates a tidy-format table with columns Timepoint_h, Organ, Passive_%ID_g, Active_%ID_g, and Enhancement, enabling direct plotting or statistical analysis in any external tool.

`export_pk_parameters_to_csv()` provides a wide-format comparison table for five formulations with all primary pharmacokinetic parameters plus an AUC fold-improvement column computed relative to free doxorubicin.

`export_toxicity_to_csv()` exports mean and standard deviation values for all toxicity biomarkers across four treatment groups.

---

## Extension Guide

### Adding a New Mathematical Model

Add the function to `src/models.py` following the docstring template used throughout the module. The docstring must include a one-line summary, the governing equation in ASCII notation, a description of each parameter with its units, the return value description, and at least one primary literature citation. Then add a corresponding unit test class in `tests/test_models.py` covering at least a shape or dimensional check, a known analytical solution, biological plausibility bounds, and relevant edge cases.

### Connecting to Experimental Data

The framework is designed to accept external datasets for parameter fitting. Pass measured concentration-time data to `scipy.optimize.curve_fit` with `two_compartment_pk` as the model function to obtain fitted pharmacokinetic parameters. Pass measured drug release data to `korsmeyer_peppas_release` to fit k and n. Pass measured dose-response data to `fit_ic50`. All model functions accept standard NumPy arrays, making integration with Pandas DataFrames straightforward through the `.to_numpy()` method.

---

*For questions on the code architecture or to report issues: github.com/Om-Physics/Nanoparticles/issues*
*Contact: om.physics7@gmail.com*
