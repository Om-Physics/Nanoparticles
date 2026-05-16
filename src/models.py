#!/usr/bin/env python3
"""
models.py
=========
Core mathematical and pharmacokinetic models used across all computational
analyses in the gold nanoparticle drug delivery research framework.

Author  : Om Jha, St. Xavier's College Kathmandu
Contact : om.physics7@gmail.com
Year    : 2026

All equations are derived from peer-reviewed literature and validated
against published experimental ranges for gold nanoparticle systems.
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit


# ══════════════════════════════════════════════════════════════════════════════
# PHARMACOKINETIC MODELS
# ══════════════════════════════════════════════════════════════════════════════

def two_compartment_pk(t, A, alpha, B, beta):
    """
    Two-compartment pharmacokinetic model for blood concentration.

    Equation
    --------
    C(t) = A * exp(-alpha * t) + B * exp(-beta * t)

    where alpha is the distribution rate constant and beta is the
    elimination rate constant.  The parameters A and B are intercept
    coefficients related to volume of distribution and clearance.

    Parameters
    ----------
    t     : array_like  Time points (hours).
    A     : float       Intercept of the distribution phase.
    alpha : float       Distribution rate constant (1/h).
    B     : float       Intercept of the elimination phase.
    beta  : float       Elimination rate constant (1/h).

    Returns
    -------
    ndarray  Blood concentration as percent injected dose.

    Literature
    ----------
    Gibaldi & Perrier, Pharmacokinetics, 2nd ed., 1982.
    Perrault et al., Nano Lett. 2009, 9, 1909-1915.
    """
    return A * np.exp(-alpha * t) + B * np.exp(-beta * t)


def epr_tumor_accumulation(t, k_in, k_out, plasma_params):
    """
    EPR-mediated tumour accumulation via first-order differential equation.

    Equation
    --------
    dC_tumour/dt = k_in * C_plasma(t) - k_out * C_tumour

    The EPR influx dominates over efflux during initial accumulation
    (k_in >> k_out), while clearance governs the late phase.

    Parameters
    ----------
    t             : array_like  Time vector (hours).
    k_in          : float       EPR influx rate constant (1/h).
    k_out         : float       Tumour efflux rate constant (1/h).
    plasma_params : tuple       (A, alpha, B, beta) for two_compartment_pk.

    Returns
    -------
    ndarray  Tumour concentration (%ID/g).

    Literature
    ----------
    Cabral et al., Nature Nanotech. 2011, 6, 815-823.
    Fang et al., Adv. Drug Deliv. Rev. 2011, 63, 136-151.
    """
    A, alpha, B, beta = plasma_params

    def ode(C_t, t_):
        C_p = two_compartment_pk(t_, A, alpha, B, beta)
        return k_in * C_p - k_out * C_t

    return odeint(ode, 0.0, t).flatten()


def pbpk_organ(t, Q, V, P, plasma_fn):
    """
    Simplified physiologically-based pharmacokinetic model for organ distribution.

    Equation
    --------
    dC_organ/dt = (Q / V) * (C_blood - C_organ / P)

    Parameters
    ----------
    t        : array_like  Time vector (hours).
    Q        : float       Organ blood flow (mL/h).
    V        : float       Organ volume (mL).
    P        : float       Organ-plasma partition coefficient.
    plasma_fn: callable    Function returning plasma concentration at time t.

    Returns
    -------
    ndarray  Organ concentration (%ID/g).

    Literature
    ----------
    Rowland et al., J. Pharmacokinet. Pharmacodyn. 2004, 31, 369-373.
    """
    def ode(C_org, t_):
        return (Q / V) * (plasma_fn(t_) - C_org / P)

    return odeint(ode, 0.0, t).flatten()


def calculate_pk_parameters(t, C):
    """
    Derive key pharmacokinetic parameters from a concentration-time profile.

    Computes Cmax, Tmax, AUC (trapezoidal), and apparent half-life from
    the terminal elimination phase.

    Parameters
    ----------
    t : array_like  Time points (hours).
    C : array_like  Corresponding concentrations (%ID/mL).

    Returns
    -------
    dict with keys: Cmax, Tmax, AUC, t_half
    """
    t = np.asarray(t, dtype=float)
    C = np.asarray(C, dtype=float)

    cmax  = float(np.max(C))
    tmax  = float(t[np.argmax(C)])
    auc   = float(np.trapz(C, t))

    # Estimate terminal half-life from last three points
    terminal_t = t[-3:]
    terminal_C = np.log(np.clip(C[-3:], 1e-12, None))
    if len(terminal_t) >= 2 and np.std(terminal_t) > 0:
        slope, _ = np.polyfit(terminal_t, terminal_C, 1)
        t_half = float(-np.log(2) / slope) if slope < 0 else np.nan
    else:
        t_half = np.nan

    return {"Cmax": cmax, "Tmax": tmax, "AUC": auc, "t_half": t_half}


# ══════════════════════════════════════════════════════════════════════════════
# DRUG RELEASE MODELS
# ══════════════════════════════════════════════════════════════════════════════

def first_order_release(t, k, M_total=100.0):
    """
    First-order drug release kinetics.

    Equation
    --------
    M(t) = M_total * (1 - exp(-k * t))

    Parameters
    ----------
    t       : array_like  Time (hours).
    k       : float       Release rate constant (1/h).
    M_total : float       Total drug loaded (%, default 100).

    Returns
    -------
    ndarray  Cumulative drug released (%).

    Literature
    ----------
    Equation standard for nanoparticle release profiling.
    """
    return M_total * (1.0 - np.exp(-k * t))


def korsmeyer_peppas_release(t, k, n, M_total=100.0):
    """
    Korsmeyer-Peppas model for drug release from nanoparticle matrices.

    Equation
    --------
    M_t / M_inf = k * t^n

    Mechanistic interpretation of n
    --------------------------------
    n < 0.45   Fickian diffusion
    0.45-0.89  Anomalous (non-Fickian) transport
    n > 0.89   Case-II (relaxation-controlled) transport

    Parameters
    ----------
    t       : array_like  Time (hours).
    k       : float       Release rate constant.
    n       : float       Release exponent.
    M_total : float       Maximum cumulative release (%).

    Returns
    -------
    ndarray  Cumulative drug released (%), clipped at M_total.

    Literature
    ----------
    Korsmeyer et al., Int. J. Pharm. 1983, 15, 25-35.
    """
    return np.clip(M_total * k * (t ** n), 0.0, M_total)


def ph_dependent_release(t, pH, k_max=0.05, pH_ref=6.0, hill=3.0, M_total=100.0):
    """
    pH-dependent first-order release with Hill-type pH sensitivity.

    Equation
    --------
    k(pH) = k_max / (1 + (pH / pH_ref)^hill)
    M(t)  = M_total * (1 - exp(-k(pH) * t))

    Captures enhanced release in acidic endosomal compartments
    (pH 5.0-5.5) versus physiological blood pH (7.4).

    Parameters
    ----------
    t       : array_like  Time (hours).
    pH      : float       Environmental pH.
    k_max   : float       Maximum release rate constant (1/h).
    pH_ref  : float       Reference pH (inflection point).
    hill    : float       Hill coefficient (cooperativity).
    M_total : float       Total drug loaded (%).

    Returns
    -------
    ndarray  Cumulative drug released (%).
    """
    k_ph = k_max / (1.0 + (pH / pH_ref) ** hill)
    return M_total * (1.0 - np.exp(-k_ph * t))


# ══════════════════════════════════════════════════════════════════════════════
# DOSE-RESPONSE MODELS
# ══════════════════════════════════════════════════════════════════════════════

def four_parameter_logistic(dose, bottom, top, IC50, hill):
    """
    Four-parameter logistic (4PL) dose-response model.

    Equation
    --------
    Response = bottom + (top - bottom) / (1 + (dose / IC50)^hill)

    This is the gold-standard model for IC50 determination in
    cytotoxicity assays.  Fitting is performed via nonlinear
    least-squares optimisation.

    Parameters
    ----------
    dose   : array_like  Dose or concentration values.
    bottom : float       Minimum response (typically 0 %).
    top    : float       Maximum response (typically 100 %).
    IC50   : float       Half-maximal inhibitory concentration.
    hill   : float       Hill slope (steepness of transition).

    Returns
    -------
    ndarray  Predicted cell viability (%).

    Literature
    ----------
    DeLean et al., Am. J. Physiol. 1978, 235, E97-E102.
    """
    return bottom + (top - bottom) / (1.0 + (dose / IC50) ** hill)


def fit_ic50(doses, responses, p0=None):
    """
    Fit the 4PL model to dose-response data and return IC50.

    Parameters
    ----------
    doses     : array_like  Dose values.
    responses : array_like  Measured cell viability (%).
    p0        : list        Initial guess [bottom, top, IC50, hill].

    Returns
    -------
    ic50    : float  Fitted IC50 value.
    popt    : ndarray  All fitted parameters [bottom, top, IC50, hill].
    pcov    : ndarray  Parameter covariance matrix.
    r2      : float  Coefficient of determination.
    """
    doses     = np.asarray(doses,     dtype=float)
    responses = np.asarray(responses, dtype=float)

    if p0 is None:
        p0 = [0.0, 100.0, np.median(doses), -1.5]

    try:
        popt, pcov = curve_fit(
            four_parameter_logistic, doses, responses, p0=p0,
            bounds=([0, 50, 1e-3, -10], [20, 100, 1e6, 0]),
            maxfev=10000
        )
    except RuntimeError:
        return np.nan, np.array(p0), np.eye(4) * np.inf, 0.0

    residuals = responses - four_parameter_logistic(doses, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((responses - np.mean(responses)) ** 2)
    r2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return float(popt[2]), popt, pcov, float(r2)


# ══════════════════════════════════════════════════════════════════════════════
# PHOTOTHERMAL MODELS
# ══════════════════════════════════════════════════════════════════════════════

def photothermal_heating(t, eta, I, tau, T_body=37.0):
    """
    Temperature rise during laser irradiation of gold nanoparticles.

    Equation
    --------
    T(t) = T_body + (eta * I) * tau * (1 - exp(-t / tau))

    The steady-state temperature rise equals eta * I * tau, reached
    asymptotically with time constant tau.

    Parameters
    ----------
    t      : array_like  Time (seconds).
    eta    : float       Photothermal conversion efficiency (0-1).
    I      : float       Laser power density (W/cm²).
    tau    : float       Thermal time constant (seconds).
    T_body : float       Initial tissue temperature (°C), default 37.

    Returns
    -------
    ndarray  Absolute temperature (°C).

    Literature
    ----------
    Cole et al., J. Phys. Chem. C 2009, 113, 12090-12094.
    """
    return T_body + (eta * I) * tau * (1.0 - np.exp(-t / tau))


def photothermal_cooling(t, T_peak, tau_cool, T_body=37.0):
    """
    Exponential tissue cooling after laser is switched off.

    Equation
    --------
    T(t) = T_body + (T_peak - T_body) * exp(-t / tau_cool)

    Parameters
    ----------
    t        : array_like  Time after laser-off (seconds).
    T_peak   : float       Peak temperature at laser-off (°C).
    tau_cool : float       Cooling time constant (seconds).
    T_body   : float       Baseline tissue temperature (°C).

    Returns
    -------
    ndarray  Absolute temperature (°C).
    """
    return T_body + (T_peak - T_body) * np.exp(-t / tau_cool)


def cem43_dose(T_profile, dt=1.0):
    """
    Cumulative Equivalent Minutes at 43 °C (CEM43) thermal dose.

    Equation
    --------
    CEM43 = sum( R^(43 - T) * dt )

    where R = 0.25 for T < 43 °C and R = 0.5 for T >= 43 °C.

    CEM43 > 60 is generally accepted as the lethal threshold for
    irreversible tissue damage.

    Parameters
    ----------
    T_profile : array_like  Temperature profile (°C).
    dt        : float       Time step (minutes).

    Returns
    -------
    ndarray  Cumulative CEM43 at each time step.

    Literature
    ----------
    Sapareto & Dewey, Int. J. Radiat. Oncol. Biol. Phys. 1984, 10, 787-800.
    """
    T = np.asarray(T_profile, dtype=float)
    R = np.where(T < 43.0, 0.25, 0.5)
    return np.cumsum(R ** (43.0 - T) * dt)


# ══════════════════════════════════════════════════════════════════════════════
# NANOPARTICLE SIZE AND OPTICAL PROPERTIES
# ══════════════════════════════════════════════════════════════════════════════

def mie_extinction_coefficient(R, epsilon_m, epsilon_r, epsilon_i, wavelength):
    """
    Mie theory extinction coefficient for spherical gold nanoparticles.

    Equation (quasi-static limit)
    --------
    epsilon = (24 * pi^2 * N_A * R^3 * epsilon_m^(3/2)) /
              (lambda * ln(10)) *
              (epsilon_i / ((epsilon_r + 2*epsilon_m)^2 + epsilon_i^2))

    Parameters
    ----------
    R          : float  Particle radius (nm).
    epsilon_m  : float  Dielectric constant of surrounding medium.
    epsilon_r  : float  Real part of gold dielectric function.
    epsilon_i  : float  Imaginary part of gold dielectric function.
    wavelength : float  Wavelength (nm).

    Returns
    -------
    float  Molar extinction coefficient (M⁻¹ cm⁻¹).

    Literature
    ----------
    Haiss et al., Anal. Chem. 2007, 79, 4215-4221.
    Link & El-Sayed, J. Phys. Chem. B 1999, 103, 4212-4217.
    """
    NA   = 6.022e23
    R_cm = R * 1e-7      # nm to cm
    lam  = wavelength * 1e-7  # nm to cm

    numerator   = 24 * np.pi**2 * NA * R_cm**3 * epsilon_m**(3/2)
    denominator = lam * np.log(10)
    lorentzian  = epsilon_i / ((epsilon_r + 2*epsilon_m)**2 + epsilon_i**2)

    return (numerator / denominator) * lorentzian


def stokes_einstein_radius(D, T=298.15, eta=0.001):
    """
    Hydrodynamic radius from Stokes-Einstein diffusion coefficient.

    Equation
    --------
    R_h = k_B * T / (6 * pi * eta * D)

    Parameters
    ----------
    D   : float  Diffusion coefficient (m²/s).
    T   : float  Temperature (K), default 298.15 (25 °C).
    eta : float  Solvent viscosity (Pa·s), default water at 25 °C.

    Returns
    -------
    float  Hydrodynamic radius (nm).

    Literature
    ----------
    Einstein, Ann. Phys. 1905, 17, 549-560.
    """
    k_B = 1.380649e-23   # J/K
    R_h = k_B * T / (6.0 * np.pi * eta * D)
    return R_h * 1e9     # m to nm


def scherrer_crystallite_size(beta_rad, theta_rad, K=0.9, wavelength_nm=0.154):
    """
    Scherrer equation for XRD-based crystallite size.

    Equation
    --------
    D = K * lambda / (beta * cos(theta))

    Parameters
    ----------
    beta_rad     : float  Full-width at half-maximum (radians).
    theta_rad    : float  Bragg angle (radians).
    K            : float  Shape factor (0.9 for spheres).
    wavelength_nm: float  X-ray wavelength in nm (Cu Kα = 0.154 nm).

    Returns
    -------
    float  Crystallite size (nm).

    Literature
    ----------
    Scherrer, Göttinger Nachrichten Gesell. 1918, 2, 98-100.
    """
    return K * wavelength_nm / (beta_rad * np.cos(theta_rad))


# ══════════════════════════════════════════════════════════════════════════════
# RECEPTOR-LIGAND BINDING
# ══════════════════════════════════════════════════════════════════════════════

def receptor_mediated_uptake(t, N_NP, N_R, k_on, k_off, k_int):
    """
    Kinetic model of receptor-mediated nanoparticle endocytosis.

    Equations
    ---------
    dN_bound/dt = k_on * N_NP * N_R - k_off * N_bound - k_int * N_bound
    dN_cell/dt  = k_int * N_bound

    Parameters
    ----------
    t     : array_like  Time (hours).
    N_NP  : float       Extracellular nanoparticle concentration.
    N_R   : float       Surface receptor density (per cell).
    k_on  : float       Association rate constant (1/(nM·h)).
    k_off : float       Dissociation rate constant (1/h).
    k_int : float       Internalisation rate constant (1/h).

    Returns
    -------
    tuple  (N_bound(t), N_cell(t)) — bound and internalised counts.

    Literature
    ----------
    Zhang et al., Adv. Mater. 2009, 21, 419-424.
    """
    def system(y, t_):
        N_bound, N_cell = y
        dbound = k_on * N_NP * N_R - (k_off + k_int) * N_bound
        dcell  = k_int * N_bound
        return [dbound, dcell]

    sol = odeint(system, [0.0, 0.0], t)
    return sol[:, 0], sol[:, 1]


def langmuir_adsorption(C, C_max, K_d):
    """
    Langmuir adsorption isotherm for drug loading on nanoparticle surface.

    Equation
    --------
    q(C) = C_max * C / (K_d + C)

    Parameters
    ----------
    C     : array_like  Free drug concentration (µg/mL).
    C_max : float       Maximum surface loading capacity (µg drug / mg Au).
    K_d   : float       Dissociation constant (µg/mL).

    Returns
    -------
    ndarray  Surface-adsorbed drug amount (µg drug / mg Au).

    Literature
    ----------
    Langmuir, J. Am. Chem. Soc. 1918, 40, 1361-1403.
    """
    return C_max * C / (K_d + C)


# ══════════════════════════════════════════════════════════════════════════════
# EPR EFFECT QUANTIFICATION
# ══════════════════════════════════════════════════════════════════════════════

def epr_accumulation_efficiency(diameter_nm, zeta_mV, t_circ_h, K_perm):
    """
    Simplified EPR accumulation efficiency model.

    Equation
    --------
    E = f(d, zeta, t_circ, K_perm)

    A heuristic combining size-dependent extravasation, surface-charge
    stability, circulation time, and vascular permeability coefficient.

    Parameters
    ----------
    diameter_nm : float  Nanoparticle hydrodynamic diameter (nm).
    zeta_mV     : float  Zeta potential (mV).
    t_circ_h    : float  Blood circulation half-life (h).
    K_perm      : float  Vascular permeability coefficient (cm/s).

    Returns
    -------
    float  Relative accumulation efficiency (0-1 scale).

    Literature
    ----------
    Fang et al., Adv. Drug Deliv. Rev. 2011, 63, 136-151.
    """
    size_factor  = np.exp(-((diameter_nm - 50) ** 2) / (2 * 30 ** 2))
    charge_factor = 1.0 - np.exp(-abs(zeta_mV) / 30.0)
    circ_factor   = 1.0 - np.exp(-t_circ_h / 20.0)
    perm_factor   = K_perm / (K_perm + 1e-7)

    return size_factor * charge_factor * circ_factor * perm_factor


# ══════════════════════════════════════════════════════════════════════════════
# STARLING EQUATION (FLUID FLUX)
# ══════════════════════════════════════════════════════════════════════════════

def starling_fluid_flux(L_p, S, P_c, P_i, sigma, pi_c, pi_i):
    """
    Starling equation for fluid flux across capillary walls.

    Equation
    --------
    J_v = L_p * S * [(P_c - P_i) - sigma * (pi_c - pi_i)]

    Describes the hydrostatic and oncotic pressure balance governing
    vascular permeability and nanoparticle extravasation.

    Parameters
    ----------
    L_p   : float  Hydraulic conductivity (cm / (mmHg·s)).
    S     : float  Capillary surface area (cm²).
    P_c   : float  Capillary hydrostatic pressure (mmHg).
    P_i   : float  Interstitial hydrostatic pressure (mmHg).
    sigma : float  Reflection coefficient (0-1).
    pi_c  : float  Capillary oncotic pressure (mmHg).
    pi_i  : float  Interstitial oncotic pressure (mmHg).

    Returns
    -------
    float  Net fluid flux J_v (mL/s).

    Literature
    ----------
    Stylianopoulos et al., PNAS 2012, 109, 15101-15108.
    """
    return L_p * S * ((P_c - P_i) - sigma * (pi_c - pi_i))
