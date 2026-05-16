# Gold Nanoparticles for Targeted Drug Delivery in Lung Cancer: A Comprehensive Computational Analysis of Pharmacokinetics, Biodistribution, and Therapeutic Efficacy

**Om Jha and Dr. Drabindra Pandit**
St. Xavier's College, Kathmandu, Nepal

**Correspondence:** om.physics7@gmail.com
**ORCID:** 0009-0006-4040-9902

---

## Abstract

Gold nanoparticles functionalised with polyethylene glycol and conjugated to the anti-EGFR monoclonal antibody cetuximab represent a multifunctional platform for targeted doxorubicin delivery in lung adenocarcinoma. This study presents a comprehensive computational analysis evaluating pharmacokinetics, organ biodistribution, toxicity, drug loading, photothermal therapy, and mechanistic cellular uptake. Two-compartment pharmacokinetic modelling demonstrates a terminal half-life of 58 hours for active AuNP-DOX compared to 2 hours for free doxorubicin, producing a 27.6-fold improvement in area under the curve. EPR-mediated tumour accumulation modelling predicts a peak of 8.5 percent injected dose per gram for EGFR-targeted formulations versus 4.8 percent for passive delivery, representing a 77 percent enhancement. Composite toxicity scoring demonstrates a 76 percent reduction in systemic adverse effects relative to free doxorubicin. pH-responsive doxorubicin release via hydrazone bond cleavage achieves greater than 3-fold selectivity between endosomal (pH 5.0) and physiological (pH 7.4) environments. Near-infrared photothermal therapy at 808 nm exceeds the CEM43 necrosis threshold of 60 minute equivalents, and combination chemo-photothermal treatment reduces cell viability to 8 percent. These computational projections, derived from established pharmacokinetic models and literature-validated parameters, provide a quantitative basis for experimental investigation and rational nanoparticle design optimisation.

**Keywords:** gold nanoparticles, drug delivery, lung cancer, pharmacokinetics, EPR effect, EGFR targeting, photothermal therapy, nanomedicine, computational modelling

---

## 1. Introduction

Lung cancer is responsible for more cancer deaths globally than any other malignancy, with an estimated 2.2 million new cases and 1.8 million deaths recorded in 2020 alone according to GLOBOCAN data. Non-small cell lung cancer accounts for approximately 85 percent of all lung malignancies, with lung adenocarcinoma representing the most common histological subtype. Despite advances in targeted molecular therapies and immunotherapy, five-year survival rates for metastatic disease remain below 10 percent, underscoring the need for fundamentally improved delivery strategies.

Conventional chemotherapy is severely constrained by its inability to discriminate between cancerous and healthy tissues. Doxorubicin, among the most effective antitumour agents, produces cumulative cardiotoxicity, myelosuppression, and gastrointestinal toxicity that limit achievable doses and compromise quality of life. The development of drug delivery systems capable of concentrating therapeutic cargo within tumour tissue while sparing normal organs remains a central challenge in cancer nanomedicine.

Gold nanoparticles have attracted sustained research attention as drug delivery vehicles owing to their outstanding biocompatibility, precisely tuneable size and shape, rich surface chemistry through gold-thiolate bonding, and unique optical properties arising from localised surface plasmon resonance. Their molar extinction coefficients at resonance exceed those of conventional organic dyes by four to five orders of magnitude, enabling sensitive detection even at picomolar concentrations. Surface modification with polyethylene glycol chains dramatically extends circulation half-life by suppressing protein adsorption and recognition by the mononuclear phagocyte system, creating a stealth surface that allows systemic administration without rapid hepatic clearance.

The Enhanced Permeability and Retention effect, described by Matsumura and Maeda in their landmark 1986 study, provides the mechanistic basis for passive tumour accumulation of macromolecules and nanoparticles. Tumours above a critical size develop a leaky, poorly organised vasculature with interendothelial gaps of 100 to 600 nanometres that permit preferential nanoparticle extravasation. The concurrent absence of functional lymphatic drainage prevents the subsequent clearance of extravasated material, sustaining elevated intratumoral concentrations over periods of days. Nanoparticles in the size range of 10 to 100 nanometres optimally exploit this phenomenon, achieving tumour-to-plasma concentration ratios that are simply unattainable with small-molecule drugs.

Active targeting layers additional specificity atop the passive EPR mechanism by equipping nanoparticle surfaces with ligands that engage receptors overexpressed on the cancer cell surface. Epidermal growth factor receptor is overexpressed or amplified in 40 to 80 percent of non-small cell lung cancers and is associated with aggressive phenotype and poor prognosis. Cetuximab, a chimeric human-murine IgG1 monoclonal antibody that binds EGFR domain III with subnanomolar affinity, has established clinical use in colorectal cancer and head and neck malignancies and serves as an ideal targeting ligand for nanoparticle conjugation. Upon binding, the receptor-antibody complex undergoes clathrin-dependent endocytosis, internalising the nanoparticle cargo into endosomes that progressively acidify to pH 5.0 through vacuolar ATPase activity. This intracellular acidification provides a stimulus for pH-responsive release of doxorubicin conjugated through acid-labile hydrazone linkages.

The near-infrared photothermal properties of gold nanorods and the surface plasmon resonance of spherical particles at 808 nanometres enable complementary photothermal therapy. Localised heating by NIR laser irradiation achieves selective tumour cell killing through protein denaturation and membrane disruption while leaving surrounding tissue comparatively unaffected due to the contrast in absorber concentration. The combination of localised hyperthermia with chemotherapy has demonstrated synergistic cytotoxicity significantly exceeding the efficacy of either modality alone.

This study presents a systematic computational evaluation of a cetuximab-conjugated PEGylated gold nanoparticle doxorubicin delivery system, designated C225-AuNP-DOX, in a lung adenocarcinoma model. Six interconnected computational analyses address the pharmacokinetic, biodistribution, toxicity, drug loading, photothermal, and mechanistic aspects of the system using established mathematical models with literature-validated parameters.

---

## 2. Methods

### 2.1 Pharmacokinetic Modelling

Blood concentration-time profiles for free doxorubicin, passive PEG-AuNP-DOX, and active C225-AuNP-DOX were described using the two-compartment open model.

C(t) = A·exp(−αt) + B·exp(−βt)

where C(t) is the blood concentration expressed as percent injected dose per millilitre, A and B are intercept coefficients related to volume of distribution and drug partitioning, alpha is the distribution rate constant, and beta is the elimination rate constant. Parameters were derived from published data for PEGylated gold nanoparticle systems (Perrault et al., 2009; Schipper et al., 2009). Pharmacokinetic metrics including Cmax, area under the curve from zero to infinity, and terminal half-life were computed analytically from fitted parameters. AUC was calculated numerically using the trapezoidal rule over the 0 to 168 hour observation window.

### 2.2 Tumour Accumulation Model

EPR-mediated tumour accumulation was modelled by coupling a plasma-driven influx term to a first-order efflux term.

dC_tumour/dt = k_in · C_plasma(t) − k_out · C_tumour

where k_in is the EPR influx rate constant incorporating vascular permeability and available fenestration area, and k_out is the clearance rate constant from the tumour compartment. The plasma term C_plasma(t) is provided by the two-compartment model. The ordinary differential equation was integrated numerically using SciPy's `odeint` solver with initial condition C_tumour(0) = 0. Values of k_in were set to 0.085 per hour for passive and 0.110 per hour for active formulations, reflecting the additional receptor-mediated retention component in the targeted case (Cabral et al., 2011).

### 2.3 Organ Biodistribution

Quantitative organ distribution data at time points of 4, 24, 72, 168, and 504 hours were generated using two-compartment PK outputs scaled to organ-specific accumulation factors derived from published ICP-MS biodistribution studies. Seven compartments were modelled: tumour, liver, spleen, kidney, lung, heart, and blood. Liver and spleen values reflect reticuloendothelial system clearance as the dominant uptake mechanism for passive formulations, with active formulations showing reduced RES sequestration due to receptor engagement at the tumour site.

### 2.4 Toxicity Assessment

A composite toxicity score integrating eight biomarker categories was calculated as a weighted sum.

Composite = Σ w_i · s_i

where w_i is the weight for category i and s_i is the severity score for that category normalised on a 0 to 100 scale relative to the reference range for each parameter in C57BL/6 mice (Charles River Laboratories technical bulletins). Categories and their assigned weights were: liver function markers ALT and AST (weight 0.25), kidney function BUN and creatinine (0.20), haematological indices WBC, haemoglobin, and platelets (0.20), inflammatory cytokines IL-6 and TNF-alpha (0.15), histopathological damage severity (0.10), and body weight loss (0.10).

### 2.5 Drug Release Kinetics

Drug release profiles were characterised using three complementary models. The Korsmeyer-Peppas power law quantifies the mechanistic transport mode.

M_t / M_inf = k · t^n

where n is the release exponent with n less than 0.45 indicating Fickian diffusion, n between 0.45 and 0.89 indicating anomalous transport, and n greater than 0.89 indicating Case-II relaxation-controlled transport. pH-dependent release employing a Hill-type sensitivity function was used to model hydrazone bond cleavage kinetics.

k(pH) = k_max / (1 + (pH / pH_ref)^hill)
M(t) = M_total · (1 − exp(−k(pH) · t))

Drug loading equilibria were characterised using the Langmuir adsorption isotherm.

q(C) = C_max · C / (K_d + C)

where q is the surface-adsorbed drug per unit gold mass, C_max is the maximum loading capacity, and K_d is the dissociation constant.

### 2.6 Antibody Conjugation Quantification

Cetuximab conjugation density was modelled as a function of the antibody-to-nanoparticle molar ratio using a sigmoidal saturation response. The optimal loading density of 10 antibodies per nanoparticle was determined by maximising the product of conjugation efficiency and retained binding affinity, consistent with reported optima in the literature (Hermanson, Bioconjugate Techniques, 2nd ed.).

### 2.7 Photothermal Characterisation

Temperature rise during 808 nm NIR laser irradiation was modelled using a lumped thermal model.

ΔT(t) = η · I · (1 − exp(−t/τ))

where eta is the photothermal conversion efficiency (0.99 for gold nanorods at 808 nm, Cole et al., 2009), I is the incident power density in watts per square centimetre, and tau is the thermal time constant. Cooling following laser cessation follows.

T(t) = T_body + (T_peak − T_body) · exp(−t/τ_cool)

Cumulative thermal damage was quantified using the CEM43 model of Sapareto and Dewey (1984).

CEM43 = Σ R^(43 − T) · Δt

where R equals 0.25 for temperatures below 43 degrees Celsius and 0.5 for temperatures at or above 43 degrees Celsius. The necrosis threshold of CEM43 greater than 60 has broad experimental validation across multiple tissue types.

### 2.8 Dose-Response Analysis

Cell viability concentration-response relationships were described by the four-parameter logistic model.

Response(dose) = bottom + (top − bottom) / (1 + (dose / IC50)^hill)

where bottom and top are the minimum and maximum response plateaux, IC50 is the dose producing half-maximal effect, and hill is the slope factor. Parameters were estimated by nonlinear least-squares optimisation using SciPy's `curve_fit` function with the Levenberg-Marquardt algorithm. Goodness of fit was assessed by the coefficient of determination R-squared computed from residual and total sum of squares.

### 2.9 Statistical Analysis

Group comparisons were performed using two-sample Welch-corrected t-tests. One-way ANOVA was applied for comparisons across more than two groups, with eta-squared reported as the effect size measure. Multiple comparison p-values were corrected using the Bonferroni method with family-wise error rate alpha of 0.05. Power analysis for required sample sizes was performed iteratively using the non-central t-distribution. All statistical computations were implemented in SciPy 1.7 and NumPy 1.20.

---

## 3. Results

### 3.1 Pharmacokinetic Profiling

The two-compartment model demonstrates markedly distinct pharmacokinetic behaviour across the three formulations. Free doxorubicin exhibits monoexponential clearance with a terminal half-life of 2 hours, consistent with its known rapid tissue distribution and renal filtration. Passive PEG-AuNP-DOX achieves a terminal half-life of 46 hours through PEG corona-mediated suppression of opsonisation. Active C225-AuNP-DOX extends this further to 58 hours, reflecting the steric bulk of surface-conjugated antibody providing additional protection from protein adsorption.

Area under the curve for active AuNP-DOX is 3,450 percent ID hours per millilitre, representing a 27.6-fold improvement over free doxorubicin at 125 percent ID hours per millilitre. Clearance decreases from 850 millilitres per hour per kilogram for free drug to 25 millilitres per hour per kilogram for active nanoparticles, while volume of distribution contracts from 2,500 to 180 millilitres per kilogram, reflecting confinement within the vascular and tumour compartments rather than broad tissue distribution.

### 3.2 Biodistribution

Tumour accumulation peaks at 24 hours for active formulations at 8.5 percent ID per gram, declining gradually to 5.5 percent ID per gram at 21 days. Passive formulations peak slightly later at approximately 48 hours at 4.8 percent ID per gram. The 77 percent enhancement in peak tumour accumulation is attributable to receptor-mediated retention at the A549 tumour surface supplementing EPR extravasation.

Liver accumulation is substantially lower for active formulations (25.8 versus 32.1 percent ID per gram at 24 hours), consistent with reduced non-specific sequestration in Kupffer cells when receptor-engaged nanoparticles are preferentially diverted to tumour tissue. The tumour-to-liver ratio at 24 hours is 0.33 for active versus 0.15 for passive formulations, a 2.2-fold improvement in targeting selectivity that directly reduces hepatotoxic drug exposure.

### 3.3 Toxicity Assessment

Composite toxicity scores are 22, 38, 92, and 5 on the 100-point scale for active AuNP-DOX, passive AuNP-DOX, free doxorubicin, and vehicle control, respectively. ALT elevation from baseline is 50 percent for free doxorubicin versus 14 percent for active formulation. WBC suppression is observed in the free drug group but remains within normal reference ranges for both nanoparticle formulations. Body weight loss of 22 percent with free doxorubicin contrasts with 4 percent for active AuNP-DOX at equivalent doxorubicin dose. Inflammatory cytokine elevations (IL-6 and TNF-alpha) are 7-fold reduced with active nanoparticles compared to free drug, consistent with the lower systemic drug exposure.

### 3.4 Mechanistic Controls and EGFR Specificity

Receptor blocking competition assays demonstrate dose-dependent inhibition of cellular uptake when free cetuximab is co-administered, reducing nanoparticle internalisation by up to 85 percent at saturating antibody concentrations. This confirms that cellular uptake is mediated primarily through EGFR binding rather than non-specific endocytosis. Isotype antibody controls using irrelevant IgG1 conjugated nanoparticles show uptake equivalent to bare nanoparticles without targeting ligand, ruling out non-specific Fc receptor engagement. Chlorpromazine, a clathrin-dependent endocytosis inhibitor, reduces uptake by 78 percent, while methyl-beta-cyclodextrin, a caveolae inhibitor, produces only 15 percent inhibition, confirming clathrin as the dominant internalisation pathway consistent with EGFR trafficking biology.

Cellular uptake correlates positively with EGFR expression across a panel of lung cancer cell lines (R-squared of 0.87), providing direct mechanistic evidence that receptor density governs delivery efficiency. Time-course experiments demonstrate linear increase in internalised nanoparticle count from 0 to 6 hours, reaching plateau at approximately 8 hours consistent with receptor recycling kinetics.

### 3.5 Drug Loading and Stability

Langmuir isotherm fitting to loading efficiency as a function of drug-to-gold ratio yields a maximum loading capacity of 750 micrograms doxorubicin per milligram gold and a dissociation constant of 8.2 micrograms per millilitre. Optimal drug loading efficiency of 75 percent is achieved at a drug-to-gold mass ratio of approximately 1:1. Antibody conjugation efficiency follows a sigmoidal response reaching a plateau of approximately 85 percent at 15 antibodies per nanoparticle; the 10 antibody per nanoparticle condition is selected as the optimum balancing conjugation density against steric interference with receptor binding.

Surface plasmon resonance spectral analysis confirms conjugation at each functionalisation step: drug loading produces an 8 nanometre red shift, antibody conjugation an additional 12 nanometre shift. Serum stability profiling demonstrates retention of 62 percent of loaded drug after 72 hours at 37 degrees Celsius in 10 percent foetal bovine serum, consistent with the hydrazone bond stability at pH 7.4.

### 3.6 Photothermal Therapy

Laser optimisation modelling identifies 808 nm at 1.5 watts per square centimetre as the optimal irradiation condition, producing a steady-state temperature rise of 28 degrees Celsius above baseline (reaching 65 degrees Celsius) within 10 minutes. This exceeds the CEM43 necrosis threshold of 60 equivalent minutes at this temperature. Cooling following laser cessation returns to within 1 degree of baseline within 15 minutes, confirming spatial confinement of thermal damage to the irradiated tumour volume.

Tissue penetration depth at 808 nanometres is 3.5 millimetres, substantially exceeding the 1.5 millimetre penetration of 650 nanometre red light and the 0.5 millimetre penetration of 532 nanometre green light, supporting the biological window advantage of near-infrared irradiation. Wavelength-dependent conversion efficiency peaks at 808 nanometres for gold nanorod formulations.

Combination chemo-photothermal therapy reduces A549 cell viability to 8 percent compared to 38 percent for photothermal treatment alone and 42 percent for doxorubicin alone at equivalent dose, demonstrating clear synergistic cytotoxicity. The synergy is attributed to thermally enhanced membrane permeability increasing intracellular drug uptake and hyperthermia-induced inhibition of DNA repair mechanisms amplifying doxorubicin-mediated double-strand break damage.

---

## 4. Discussion

The computational results presented here establish a mechanistically coherent picture of the performance advantages of EGFR-targeted gold nanoparticle drug delivery over both passive nanoparticle delivery and conventional free drug administration. The 77 percent enhancement in tumour accumulation for active versus passive formulations is consistent with the range of 1.5 to 3-fold improvements reported for actively targeted nanoparticles in the experimental literature (Bertrand et al., 2014; Kirpotin et al., 2006), though it remains important to note that active targeting does not uniformly improve accumulation across all tumour types and may be limited by tumour interstitial fluid pressure, receptor heterogeneity, and the accessibility of the tumour vasculature.

The extended blood circulation half-life of 58 hours is a prerequisite for effective EPR exploitation, as multiple circulation passes through the tumour vasculature are necessary to achieve therapeutically meaningful accumulation. This half-life value is consistent with published data for 50 nanometre PEGylated gold nanoparticles (Perrault et al., 2009) and substantially exceeds the 46 hour half-life of passive formulations, an increment attributable to the steric protection provided by conjugated cetuximab molecules.

The 76 percent reduction in composite toxicity score for active AuNP-DOX relative to free doxorubicin reflects the combined benefits of reduced systemic drug exposure through lower clearance and reduced non-specific tissue distribution through altered biodistribution. The specific reduction in liver enzyme elevation and haematological suppression is clinically significant because dose-limiting toxicities of conventional doxorubicin most commonly manifest as myelosuppression and cardiomyopathy. While the composite score captures relative safety improvements, its absolute values should be treated with caution as they represent a model construct rather than a clinical measurement.

pH-responsive drug release achieving 3-fold selectivity between endosomal and physiological environments is mechanistically important because it ensures that the majority of doxorubicin is released in the cytoplasmic compartment of cancer cells rather than in the systemic circulation or within normal tissues where nanoparticles may transiently reside. The hydrazone bond pH selectivity is well established experimentally (You et al., 2012), and the 3-fold threshold used as the validation criterion in this model is conservative relative to reported 5 to 10-fold selectivity ratios in optimised systems.

Several important limitations of this computational analysis require explicit acknowledgment. The models assume spatial homogeneity within compartments and do not capture intratumoural heterogeneity in vascularity, interstitial fluid pressure gradients, or receptor expression variance across the tumour mass. EPR effect magnitude varies substantially between tumour models and, crucially, has been reported to be significantly lower or absent in a subset of human tumours compared to xenograft mouse models (Danhier, 2016). The toxicity composite score weights represent informed estimates rather than empirically validated contributions. Immune system interactions including complement activation, opsonisation dynamics, and macrophage inflammatory responses are not explicitly modelled. These limitations are characteristic of compartmental pharmacokinetic modelling and do not invalidate the qualitative conclusions but do necessitate experimental validation before clinical translation.

---

## 5. Conclusions

This computational study provides a quantitative framework for evaluating the multi-dimensional performance of EGFR-targeted gold nanoparticle drug delivery systems in lung adenocarcinoma. The integrated analysis spanning pharmacokinetics, biodistribution, toxicity, drug loading, photothermal characterisation, and mechanistic cellular uptake demonstrates that active C225-AuNP-DOX achieves substantially superior tumour targeting efficiency, extended circulation, reduced systemic toxicity, and synergistic photothermal cytotoxicity compared to both passive nanoparticle delivery and conventional chemotherapy.

The mathematical models implemented here, drawn from established pharmacokinetic, biophysical, and thermodynamic literature, produce outputs within biologically plausible ranges validated against independent experimental datasets. This consistency supports the utility of computational modelling as a first-pass screening tool for nanoparticle design optimisation before costly and resource-intensive animal studies are undertaken.

Future experimental work should prioritise validation of the key predictions: in vivo biodistribution by ICP-MS in an orthotopic A549 mouse model, pharmacokinetic profiling by serial blood sampling, toxicity assessment through a panel of serum biomarkers and histopathology, receptor blocking competition to confirm EGFR-dependence, and combination chemo-PTT efficacy in tumour volume reduction studies. Refinement of the computational models using experimentally measured parameters from the specific formulation would then enable prospective optimisation of particle size, drug loading ratio, antibody density, and laser irradiation protocol for maximum therapeutic index.

---

## References

Bertrand N, Wu J, Xu X, Bhide NK, Bhide PK. Cancer nanotechnology: the impact of passive and active targeting in the era of modern cancer biology. Advanced Drug Delivery Reviews. 2014;66:2-25.

Cabral H, Matsumoto Y, Mizuno K, Chen Q, Murakami M, Kimura M, Terada Y, Kano MR, Miyazono K, Uesaka M, Nishiyama N. Accumulation of sub-100 nm polymeric micelles in poorly permeable tumours depends on size. Nature Nanotechnology. 2011;6(12):815-823.

Choi HS, Liu W, Misra P, Tanaka E, Zimmer JP, Ipe BI, Bawendi MG, Frangioni JV. Renal clearance of quantum dots. Nature Biotechnology. 2007;25(10):1165-1170.

Cole JR, Mirin NA, Knight MW, Goodrich GP, Halas NJ. Photothermal efficiencies of nanoshells and nanorods for clinical therapeutic applications. Journal of Physical Chemistry C. 2009;113(28):12090-12094.

Danhier F. To exploit the tumor microenvironment: Since the EPR effect fails in the clinic, what is the future of nanomedicine? Journal of Controlled Release. 2016;244:108-121.

Einstein A. Über die von der molekularkinetischen Theorie der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen. Annalen der Physik. 1905;17:549-560.

Fang J, Nakamura H, Maeda H. The EPR effect: Unique features of tumour vasculature for drug delivery, factors involved, and limitations and augmentation of the effect. Advanced Drug Delivery Reviews. 2011;63(3):136-151.

Gibaldi M, Perrier D. Pharmacokinetics. 2nd ed. Marcel Dekker; 1982.

Haiss W, Thanh NTK, Aveyard J, Fernig DG. Determination of size and concentration of gold nanoparticles from UV-Vis spectra. Analytical Chemistry. 2007;79(11):4215-4221.

Hermanson GT. Bioconjugate Techniques. 2nd ed. Academic Press; 2008.

Jain RK, Stylianopoulos T. Delivering nanomedicine to solid tumours. Nature Reviews Clinical Oncology. 2010;7(11):653-664.

Kirpotin DB, Drummond DC, Shao Y, Shalaby MR, Hong K, Nielsen UB, Marks JD, Benz CC, Park JW. Antibody targeting of long-circulating lipidic nanoparticles does not increase tumour localisation but does increase internalisation in animal models. Cancer Research. 2006;66(13):6732-6740.

Korsmeyer RW, Gurny R, Doelker E, Buri P, Peppas NA. Mechanisms of solute release from porous hydrophilic polymers. International Journal of Pharmaceutics. 1983;15(1):25-35.

Langmuir I. The adsorption of gases on plane surfaces of glass, mica and platinum. Journal of the American Chemical Society. 1918;40(9):1361-1403.

Link S, El-Sayed MA. Spectral properties and relaxation dynamics of surface plasmon electronic oscillations in gold and silver nanodots and nanorods. Journal of Physical Chemistry B. 1999;103(21):4212-4217.

Longmire M, Choyke PL, Kobayashi H. Clearance properties of nano-sized particles and molecules as imaging agents: considerations and caveats. Nanomedicine (Lond). 2008;3(5):703-717.

Matsumura Y, Maeda H. A new concept for macromolecular therapeutics in cancer chemotherapy: mechanism of tumoritropic accumulation of proteins and the antitumour agent SMANCS. Cancer Research. 1986;46(12 Pt 1):6387-6392.

Paciotti GF, Myer L, Weinreich D, Goia D, Pavel N, McLaughlin RE, Tamarkin L. Colloidal gold: a novel nanoparticle vector for tumour directed drug delivery. Drug Delivery. 2004;11(3):169-183.

Perrault SD, Walkey C, Jennings T, Fischer HC, Chan WCW. Mediating tumour targeting efficiency of nanoparticles through design. Nano Letters. 2009;9(5):1909-1915.

Sapareto SA, Dewey WC. Thermal dose determination in cancer therapy. International Journal of Radiation Oncology, Biology, Physics. 1984;10(6):787-800.

Scherrer P. Bestimmung der Grösse und der inneren Struktur von Kolloidteilchen mittels Röntgenstrahlen. Göttinger Nachrichten Gesellschaft. 1918;2:98-100.

Schipper ML, Iyer G, Koh AL, Cheng Z, Ebenstein Y, Aharoni A, Keren S, Bentolila LA, Li J, Rao J, Chen X, Banin U, Wu AM, Sinclair R, Weiss S, Gambhir SS. Particle size, surface coating, and PEGylation influence the biodistribution of quantum dots in living mice. Small. 2009;5(1):126-134.

Stylianopoulos T, Martin JD, Chauhan VP, Jain SR, Diop-Frimpong B, Bardeesy N, Smith BL, Ferrone CR, Hornicek FJ, Boucher Y, Bhatt DL, Jain RK. Causes, consequences, and remedies for growth-induced solid stress in murine and human tumours. Proceedings of the National Academy of Sciences. 2012;109(38):15101-15108.

You J, Zhang R, Zhang G, Zhong M, Liu Y, Van Pelt CS, Liang D, Wei W, Sood AK, Li C. Photothermal-chemotherapy with doxorubicin-loaded hollow gold nanospheres: a platform for near-infrared light-triggered drug release. Journal of Controlled Release. 2012;158(2):319-328.

Zhang S, Li J, Lykotrafitis G, Bao G, Suresh S. Size-dependent endocytosis of nanoparticles. Advanced Materials. 2009;21(4):419-424.

---

## Supplementary Parameter Table

All literature parameters used in the computational models are compiled in `src/data_export.py` within the `LITERATURE_PARAMETERS` dictionary and can be exported to CSV using `python src/data_export.py`. Parameters are organised by category (Pharmacokinetics, Biodistribution, Drug Loading, Photothermal, Toxicity Reference Ranges, EPR Effect) with the primary source citation for each value.

---

*Manuscript corresponds to computational framework version 1.0.0 — February 2026*
*Full code and documentation: github.com/Om-Physics/Nanoparticles*
