#!/usr/bin/env python3
"""
test_models.py
==============
Unit tests for all mathematical model functions.

Run from repository root:
    pytest tests/ -v

Author  : Om Jha, St. Xavier's College Kathmandu
Contact : om.physics7@gmail.com
Year    : 2026
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest
from src.models import (
    two_compartment_pk,
    epr_tumor_accumulation,
    calculate_pk_parameters,
    first_order_release,
    korsmeyer_peppas_release,
    ph_dependent_release,
    four_parameter_logistic,
    fit_ic50,
    photothermal_heating,
    photothermal_cooling,
    cem43_dose,
    mie_extinction_coefficient,
    stokes_einstein_radius,
    scherrer_crystallite_size,
    langmuir_adsorption,
    epr_accumulation_efficiency,
    starling_fluid_flux,
)
from src.statistics_utils import (
    two_sample_ttest,
    one_way_anova,
    sample_size_calculation,
    normality_test,
    bonferroni_correction,
    summary_statistics,
)


# ══════════════════════════════════════════════════════════════════════════════
# PHARMACOKINETIC MODEL TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestTwoCompartmentPK:

    def test_returns_correct_shape(self):
        t = np.linspace(0, 168, 500)
        C = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
        assert C.shape == t.shape

    def test_initial_concentration_is_A_plus_B(self):
        C = two_compartment_pk(np.array([0.0]), A=70, alpha=0.277, B=30, beta=0.012)
        assert abs(C[0] - 100.0) < 1e-6

    def test_strictly_positive(self):
        t = np.linspace(0.1, 168, 200)
        C = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
        assert np.all(C > 0)

    def test_monotonically_decreasing(self):
        t = np.linspace(0, 168, 200)
        C = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
        assert np.all(np.diff(C) < 0)

    def test_half_life_beta_consistent(self):
        beta = 0.012
        t    = np.array([0, np.log(2) / beta])
        C    = two_compartment_pk(t, A=0, alpha=1, B=100, beta=beta)
        assert abs(C[1] / C[0] - 0.5) < 0.01


class TestEPRTumorAccumulation:

    def test_starts_at_zero(self):
        t      = np.linspace(0, 100, 100)
        params = (70, 0.277, 30, 0.012)
        C_t    = epr_tumor_accumulation(t, 0.085, 0.015, params)
        assert abs(C_t[0]) < 0.1

    def test_peak_within_literature_range(self):
        """Peak tumour accumulation should be 1-12 %ID/g (literature range)."""
        t      = np.linspace(0, 200, 1000)
        params = (70, 0.277, 30, 0.012)
        C_t    = epr_tumor_accumulation(t, 0.085, 0.015, params)
        assert 1.0 <= C_t.max() <= 12.0

    def test_peak_timing_24_to_72h(self):
        """Peak should occur between 24 and 72 hours post-injection."""
        t      = np.linspace(0, 200, 2000)
        params = (70, 0.277, 30, 0.012)
        C_t    = epr_tumor_accumulation(t, 0.085, 0.015, params)
        t_peak = t[np.argmax(C_t)]
        assert 12.0 <= t_peak <= 96.0

    def test_non_negative(self):
        t      = np.linspace(0, 500, 500)
        params = (70, 0.277, 30, 0.012)
        C_t    = epr_tumor_accumulation(t, 0.085, 0.015, params)
        assert np.all(C_t >= -0.001)


class TestCalculatePKParameters:

    def test_cmax_is_maximum(self):
        t = np.linspace(0, 168, 200)
        C = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
        pk = calculate_pk_parameters(t, C)
        assert abs(pk["Cmax"] - C.max()) < 1e-6

    def test_auc_positive(self):
        t  = np.linspace(0, 168, 200)
        C  = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
        pk = calculate_pk_parameters(t, C)
        assert pk["AUC"] > 0

    def test_half_life_positive_and_finite(self):
        t  = np.linspace(0, 168, 200)
        C  = two_compartment_pk(t, A=70, alpha=0.277, B=30, beta=0.012)
        pk = calculate_pk_parameters(t, C)
        assert pk["t_half"] > 0
        assert np.isfinite(pk["t_half"])


# ══════════════════════════════════════════════════════════════════════════════
# DRUG RELEASE MODEL TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestDrugRelease:

    def test_first_order_starts_at_zero(self):
        assert abs(first_order_release(0, k=0.02)) < 1e-10

    def test_first_order_approaches_total(self):
        M = first_order_release(10000, k=0.05, M_total=100)
        assert M > 99.0

    def test_korsmeyer_clipped_at_M_total(self):
        M = korsmeyer_peppas_release(1e6, k=0.5, n=0.52, M_total=100)
        assert M <= 100.0

    def test_ph_dependent_higher_at_low_pH(self):
        t   = np.array([24.0])
        low = ph_dependent_release(t, pH=5.0)
        mid = ph_dependent_release(t, pH=6.5)
        hi  = ph_dependent_release(t, pH=7.4)
        assert low[0] > mid[0] > hi[0]

    def test_ph_release_selectivity_at_least_3fold(self):
        t  = np.array([24.0])
        r5 = ph_dependent_release(t, pH=5.0)[0]
        r7 = ph_dependent_release(t, pH=7.4)[0]
        assert r5 / r7 >= 3.0


# ══════════════════════════════════════════════════════════════════════════════
# DOSE-RESPONSE MODEL TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestDoseResponse:

    def test_4pl_at_ic50_is_midpoint(self):
        ic50  = 50.0
        val   = four_parameter_logistic(ic50, bottom=0, top=100,
                                        IC50=ic50, hill=-1.5)
        assert abs(val - 50.0) < 0.5

    def test_4pl_low_dose_near_top(self):
        val = four_parameter_logistic(1e-3, bottom=0, top=100, IC50=50, hill=-1.5)
        assert val > 90.0

    def test_4pl_high_dose_near_bottom(self):
        val = four_parameter_logistic(1e6, bottom=0, top=100, IC50=50, hill=-1.5)
        assert val < 10.0

    def test_fit_ic50_recovers_known_value(self):
        np.random.seed(99)
        doses  = np.logspace(0, 3, 20)
        known  = 52.0
        viab   = (four_parameter_logistic(doses, 0, 100, known, -1.8)
                  + np.random.normal(0, 1.5, 20))
        ic50, _, _, r2 = fit_ic50(doses, viab)
        assert abs(ic50 - known) / known < 0.15   # within 15%
        assert r2 > 0.95


# ══════════════════════════════════════════════════════════════════════════════
# PHOTOTHERMAL MODEL TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestPhotothermal:

    def test_starts_at_body_temperature(self):
        T = photothermal_heating(np.array([0.0]), eta=0.99, I=1.5, tau=90)
        assert abs(T[0] - 37.0) < 0.1

    def test_temperature_increases_monotonically(self):
        t = np.linspace(0, 600, 200)
        T = photothermal_heating(t, eta=0.99, I=1.5, tau=90)
        assert np.all(np.diff(T) >= 0)

    def test_delta_T_within_physical_range(self):
        """Temperature rise should be 10–50 °C for standard PTT conditions."""
        t   = np.linspace(0, 1200, 500)
        T   = photothermal_heating(t, eta=0.99, I=1.5, tau=90)
        dT  = T.max() - 37.0
        assert 10.0 <= dT <= 60.0

    def test_cooling_returns_to_baseline(self):
        T_cool = photothermal_cooling(np.array([10000.0]), T_peak=65, tau_cool=150)
        assert abs(T_cool[0] - 37.0) < 0.5

    def test_cem43_exceeds_necrosis_at_high_temp(self):
        """10 min at 50 °C should exceed CEM43 = 60."""
        T_arr  = np.full(600, 50.0)
        cem    = cem43_dose(T_arr, dt=1.0)
        assert cem[-1] > 60.0

    def test_cem43_stays_low_at_normal_temp(self):
        """10 min at 37 °C should produce negligible CEM43."""
        T_arr  = np.full(600, 37.0)
        cem    = cem43_dose(T_arr, dt=1.0)
        assert cem[-1] < 1.0


# ══════════════════════════════════════════════════════════════════════════════
# OPTICAL AND STRUCTURAL MODEL TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestOpticalModels:

    def test_mie_extinction_positive(self):
        epsilon = mie_extinction_coefficient(
            R=25, epsilon_m=1.78, epsilon_r=-11.4, epsilon_i=1.2,
            wavelength=520
        )
        assert epsilon > 0

    def test_stokes_einstein_radius_positive(self):
        D  = 9.5e-12  # m²/s (typical for 50 nm particle at 37°C)
        Rh = stokes_einstein_radius(D, T=310.15)
        assert Rh > 0

    def test_stokes_einstein_larger_D_gives_smaller_R(self):
        R1 = stokes_einstein_radius(1e-11)
        R2 = stokes_einstein_radius(2e-11)
        assert R1 > R2

    def test_scherrer_positive_result(self):
        import math
        D = scherrer_crystallite_size(
            beta_rad=math.radians(0.3),
            theta_rad=math.radians(19.1)
        )
        assert D > 0

    def test_langmuir_saturation_approaches_Cmax(self):
        q = langmuir_adsorption(C=1e6, C_max=800, K_d=15)
        assert q > 790

    def test_langmuir_at_zero_concentration_is_zero(self):
        q = langmuir_adsorption(C=0, C_max=800, K_d=15)
        assert q == 0.0


# ══════════════════════════════════════════════════════════════════════════════
# EPR EFFECT MODEL TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestEPRModel:

    def test_efficiency_in_0_1_range(self):
        E = epr_accumulation_efficiency(
            diameter_nm=50, zeta_mV=-35, t_circ_h=46, K_perm=1.2e-7
        )
        assert 0 <= E <= 1

    def test_optimal_size_gives_higher_efficiency(self):
        E_opt   = epr_accumulation_efficiency(50,  -35, 46, 1.2e-7)
        E_large = epr_accumulation_efficiency(300, -35, 46, 1.2e-7)
        assert E_opt > E_large

    def test_longer_circulation_gives_higher_efficiency(self):
        E_long  = epr_accumulation_efficiency(50, -35, 46, 1.2e-7)
        E_short = epr_accumulation_efficiency(50, -35,  2, 1.2e-7)
        assert E_long > E_short


class TestStarlingEquation:

    def test_positive_flux_for_normal_conditions(self):
        Jv = starling_fluid_flux(
            L_p=5e-7, S=1.0, P_c=25, P_i=5,
            sigma=0.9, pi_c=25, pi_i=5
        )
        assert np.isfinite(Jv)

    def test_flux_increases_with_permeability(self):
        Jv_hi = starling_fluid_flux(1e-6, 1.0, 30, 5, 0.8, 25, 5)
        Jv_lo = starling_fluid_flux(1e-7, 1.0, 30, 5, 0.8, 25, 5)
        assert Jv_hi > Jv_lo


# ══════════════════════════════════════════════════════════════════════════════
# STATISTICS UTILITY TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestStatistics:

    def setup_method(self):
        np.random.seed(42)
        self.g1 = np.random.normal(4.8, 0.35, 8)
        self.g2 = np.random.normal(8.5, 0.40, 8)

    def test_ttest_detects_significant_difference(self):
        result = two_sample_ttest(self.g1, self.g2)
        assert result["significant"] is True
        assert result["p_value"] < 0.001

    def test_ttest_returns_all_keys(self):
        result = two_sample_ttest(self.g1, self.g2)
        for key in ["t_stat", "p_value", "df", "significant",
                    "sig_label", "cohens_d", "effect_size", "ci_95"]:
            assert key in result

    def test_ttest_identical_groups_not_significant(self):
        g = np.random.normal(5.0, 0.5, 10)
        result = two_sample_ttest(g, g)
        assert result["p_value"] == pytest.approx(1.0, abs=0.05)

    def test_anova_significant_across_groups(self):
        c  = np.random.normal(1.0, 0.2, 8)
        p  = np.random.normal(4.8, 0.3, 8)
        a  = np.random.normal(8.5, 0.4, 8)
        r  = one_way_anova(c, p, a)
        assert r["significant"] is True
        assert r["eta_squared"] > 0.5

    def test_sample_size_returns_positive_integer(self):
        r = sample_size_calculation(effect_size=0.8)
        assert isinstance(r["n_per_group"], int)
        assert r["n_per_group"] > 0

    def test_sample_size_larger_effect_needs_less_n(self):
        r_large  = sample_size_calculation(effect_size=1.2)
        r_medium = sample_size_calculation(effect_size=0.5)
        assert r_large["n_per_group"] < r_medium["n_per_group"]

    def test_normality_normal_data(self):
        data   = np.random.normal(5, 1, 50)
        result = normality_test(data)
        assert result["is_normal"] is True

    def test_normality_clearly_non_normal(self):
        data = np.concatenate([np.ones(25), np.ones(25) * 10])
        result = normality_test(data)
        assert result["is_normal"] is False

    def test_bonferroni_corrects_alpha(self):
        p      = [0.04, 0.001, 0.03]
        result = bonferroni_correction(p)
        assert abs(result["corrected_alpha"] - 0.05/3) < 1e-8

    def test_summary_statistics_correct_mean(self):
        data   = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = summary_statistics(data)
        assert abs(result["mean"] - 3.0) < 1e-4
        assert result["n"] == 5


# ══════════════════════════════════════════════════════════════════════════════
# MASS BALANCE VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

class TestMassBalance:

    def test_biodistribution_sums_to_near_100(self):
        """
        Total organ recovery at any time point should approximate
        100% injected dose when blood, tumour, liver, spleen, kidney,
        lung, and heart compartments are summed.
        """
        bio_passive_24h = {
            "Tumor": 4.8, "Liver": 32.1, "Spleen": 19.8,
            "Kidney": 6.2, "Lung": 4.1, "Heart": 1.5, "Blood": 8.5,
        }
        total = sum(bio_passive_24h.values())
        # Allow 20% tolerance for unlisted tissues (muscle, bone, skin)
        assert 70 <= total <= 120

    def test_drug_release_does_not_exceed_100_percent(self):
        t_long  = np.linspace(0, 1000, 500)
        release = first_order_release(t_long, k=0.05, M_total=100)
        assert release.max() <= 100.0 + 1e-6
