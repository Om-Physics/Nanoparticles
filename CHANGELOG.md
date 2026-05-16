# Changelog

## Version History and Release Notes

All notable changes to this project are documented in this file following semantic versioning (MAJOR.MINOR.PATCH), where MAJOR indicates breaking changes to the API or data format, MINOR indicates new backward-compatible features or figures, and PATCH indicates bug fixes and documentation corrections.

---

### Added

**Six Computational Figures.** All six publication figures were generated, each addressing a specific reviewer concern raised during manuscript evaluation. Figure 1 provides ICP-MS-equivalent organ biodistribution data across five time points for passive and active targeting formulations, modelled using two-compartment pharmacokinetics and EPR accumulation differential equations. Figure 2 presents an eight-component safety assessment panel covering liver function, kidney function, haematology, inflammatory markers, histopathological scoring, body weight monitoring, and a composite toxicity score. Figure 3 validates EGFR targeting specificity through isotype antibody controls, dose-dependent receptor blocking competition assays, endocytosis pathway inhibitor studies, time-dependent uptake kinetics, EGFR expression correlation across cell lines, and pathway mechanism distribution. Figure 4 covers drug loading efficiency optimisation, antibody conjugation quantification, UV-visible spectral shift confirmation, serum stability profiling, BCA protein assay calibration, and fluorescence quenching analysis. Figure 5 documents laser parameter optimisation, heating and cooling cycle dynamics, wavelength-dependent photothermal conversion efficiency, chemo-PTT combination synergy, CEM43 thermal dose mapping, and tissue penetration depth comparison across five laser wavelengths. Figure 6 presents blood concentration-time curves, area under the curve comparisons with fold-improvement annotations, size-dependent clearance route analysis, and a normalised pharmacokinetic parameter summary.

**Interactive HTML5 Simulation.** A self-contained browser-based simulation was developed providing real-time visualisation of the complete eight-phase drug delivery mechanism from blood circulation through apoptosis induction. The simulation implements Brownian motion physics, receptor-ligand binding kinetics, endosomal acidification, pH-responsive drug release, nuclear translocation, and apoptosis visualisation with interactive controls, live metrics, and phase progression tracking. No external dependencies, internet connection, or server configuration are required.

**Four Core Python Modules.** The `generate_figures` module serves as the primary entry point with a command-line interface supporting selective or batch figure generation. The `models` module implements all mathematical equations including two-compartment pharmacokinetics, EPR accumulation ODEs, Korsmeyer-Peppas release kinetics, pH-dependent rate constants, Mie extinction theory, Scherrer crystallite sizing, Stokes-Einstein hydrodynamic radius, Langmuir adsorption isotherms, the Starling fluid flux equation, receptor-mediated uptake kinetics, and photothermal heating with CEM43 calculation. The `statistics_utils` module provides two-sample t-tests, one-way ANOVA, power analysis, Shapiro-Wilk normality testing, log-rank survival comparison, Bonferroni correction, and comprehensive descriptive statistics. The `data_export` module provides CSV export for all computational datasets and maintains the complete `LITERATURE_PARAMETERS` reference dictionary with source citations for every model parameter.

**Four Documentation Files.** Documentation covers the full scientific manuscript with abstract, introduction, methods, results, discussion, conclusions, and references; a complete code guide with function-level technical documentation; a figures guide with panel descriptions, reproduction instructions, and modification examples; and a simulation guide covering architecture, phase descriptions, browser compatibility, and presentation usage.

---

## Planned Future Releases

**v1.1.0 (Planned).** A Jupyter notebook tutorial series will demonstrate basic usage, parameter fitting to experimental data, sensitivity analysis, and uncertainty quantification. The biodistribution model will be extended to incorporate a full physiologically-based pharmacokinetic framework with explicit organ compartments and blood flow autoregulation. A tumour growth inhibition module will enable simulation of therapeutic response trajectories over multi-dose treatment regimens.

**v1.2.0 (Planned).** Immune system interaction modelling will be incorporated, including complement activation, macrophage recognition kinetics, and PEG anti-PEG antibody formation effects on pharmacokinetics following repeated dosing. A patient stratification module for EPR effect prediction based on tumour vascularity imaging biomarkers is planned. An optimisation framework using `scipy.optimize` will identify AuNP design parameters that maximise tumour delivery while minimising composite toxicity score.

**v2.0.0 (Long-Term).** Integration with experimental validation data from wet laboratory collaborations will enable model refinement and parameter fitting for specific synthesised formulations. A web-based interactive dashboard version of the simulation with additional visualisation modes is planned. Translation of the Python framework into a documented R package for the biostatistics community is also under consideration.

---

## Contributing Guidelines

Contributions to the codebase through bug fixes, feature additions, documentation improvements, or experimental validation data are welcome through the standard GitHub pull request process.

**Bug Reports and Fixes.** If you discover an error in any mathematical model, an incorrect parameter value, a plotting bug, or an issue with the simulation, please open a GitHub issue with sufficient detail to reproduce the problem, including your Python version, operating system, and the specific function or figure affected. For bug fixes, submit a pull request with a minimal code change targeting only the identified error, accompanied by a test case demonstrating the fix.

**New Mathematical Models.** If you wish to contribute additional pharmacokinetic, biophysical, or biological models, add the function to `src/models.py` following the established docstring format. The docstring must include a one-line summary, the governing equation in ASCII notation, parameter descriptions with units, the return value description, and at least one primary literature citation. Add a corresponding unit test in `tests/test_models.py` verifying the function returns biologically plausible values for representative inputs and correctly handles edge cases.

**New Figure Modules.** Additional figures covering topics such as tumour growth inhibition, dose-response comparison across cell lines, or receptor expression mapping are welcome. Create a self-contained function named `figure_yourtopic()` in `src/generate_figures.py`, add it to the `FIGURE_MAP` dictionary, and document it in `docs/FIGURES_GUIDE.md`. Ensure the figure uses the global `PALETTE` and `rcParams` settings for visual consistency.

**Experimental Validation Data.** Researchers who have conducted experimental studies confirming or refining the computational predictions are encouraged to contribute their data for model refinement. Contact the authors directly at om.physics7@gmail.com with a description of the experimental system, methods, and key quantitative results.

**Submission Process.** Fork the repository, create a descriptive branch name using the format `feature/description` for new features or `fix/description` for bug fixes, implement changes with clear commit messages, and submit a pull request to the main branch with a description of the contribution and relevant context. Pull requests are reviewed for scientific accuracy, code style consistency, completeness of documentation and unit tests, and correctness of literature citations. Reviews aim to be completed within two weeks.

**Code Style.** Python code follows PEP 8 conventions. Function names use lowercase with underscores. Physical quantities include their units in the variable name or docstring. Magic numbers are replaced with named constants. Functions are focused on a single mathematical or data manipulation task and do not exceed 80 lines excluding docstrings and comments. Mathematical expressions in code should closely mirror their published source equations, with comments referencing the equation number or page when applicable.

---

**Primary Contact:** om.physics7@gmail.com
**ORCID:** 0009-0006-4040-9902
**Repository:** github.com/Om-Physics/Nanoparticles
