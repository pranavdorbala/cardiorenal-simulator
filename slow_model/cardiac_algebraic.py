"""
Beat-averaged steady-state cardiac hemodynamics.

Replaces the sub-second cardiac cycle with algebraic relationships.
Given the slow state variables and parameters, computes consistent
values for all fast (beat-level) cardiac variables at quasi-steady-state.

These are *approximations* — the true beat-level dynamics involve
pulsatile flows and time-varying chamber pressures. The QSS values
are good enough for the slow model RHS evaluation and for warm-starting
the full model reconstruction.
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hallow_c_driver import PARAM_MAP, STATE_MAP


def compute_qss_fast_variables(full_y, params):
    """Compute beat-averaged values for all fast variables.

    Given a 70-dim state vector where the slow variables are meaningful
    and the fast variables may be stale, compute consistent QSS values
    for the fast variables.

    This uses the relationships from the C model equations, approximating
    all beat-level dynamics at their cycle-averaged steady state.

    Args:
        full_y: 70-dim state vector (slow vars must be valid)
        params: 430-dim parameter array

    Returns:
        qss: dict mapping state index -> QSS value for each fast variable
    """
    p = params
    y = full_y

    # ---- Myocyte geometry (from slow vars y[9], y[10]) ----
    V_w_0 = p[49]
    baseline_total_myocyte_volume = (V_w_0 - p[96] - p[97] - p[98])
    baseline_single_myocyte_volume = baseline_total_myocyte_volume / p[81]
    baseline_myocyte_diameter = 2 * np.sqrt(
        baseline_single_myocyte_volume / (p[16] * p[82]))

    myocyte_length = p[82] + y[9]
    myocyte_diameter = baseline_myocyte_diameter + y[10]
    single_myocyte_volume = myocyte_length * p[16] * myocyte_diameter**2 / 4
    total_myocyte_volume = single_myocyte_volume * p[81]
    total_nonmyocyte_volume = p[96] + p[98] + p[97]
    LV_wall_volume = total_myocyte_volume + total_nonmyocyte_volume

    # ---- Heart rate (from slow var y[50] for rsna_delayed) ----
    BB_signal = p[370] * (1 - np.exp(-p[371] * y[12]))
    beta_blocker_HR = 1 - (1 - p[368]) * BB_signal
    rsna_HR_intercept = 1 - p[298]
    rsna_effect_on_HR = p[298] * p[288] + rsna_HR_intercept
    heart_rate = p[64] * rsna_effect_on_HR * beta_blocker_HR

    # ---- LV cavity volume (depends on myocyte geometry changes) ----
    LV_cavity_volume = (p[47] *
                        (1 + p[92] * y[9] / p[82])**3 *
                        (1 - p[93] * y[10] / baseline_myocyte_diameter)**2)

    # ---- Use current CO and pressure states as QSS targets ----
    # At QSS, the delayed filters have converged, so the "old" values
    # equal the current values. We use the slow CO and existing pressures.
    CO = y[26]  # slow variable: CO in L/min
    CO_m3_per_hr = CO * p[4] / 60  # L/min -> m^3/hr (p[4] = L_m3 = 0.001)

    # Mean aortic valve flow at QSS = CO
    aortic_flow_m3_hr = CO_m3_per_hr

    # ---- Arterial compliance and pressure ----
    C_art_base = p[38]
    C_art_scale = p[116]
    # C_art transitions from C_art_base to C_art_base*C_art_scale over time
    C_art = (C_art_base - C_art_base * C_art_scale) / np.exp(
        y[12] / (24 / 2)) + C_art_base * C_art_scale

    # MAP from delayed systolic/diastolic pressures (these are slow-filtered)
    MAP_internal = y[20] / 3 + y[21] * 2 / 3  # in Pa (internal units)
    MAP_mmHg = MAP_internal * p[14]

    # Stiffness-dependent compliance
    Stiffness0 = 1 / C_art
    arterial_stiffness = Stiffness0 * (1 + (MAP_mmHg - p[126]) * p[115])
    arterial_compliance = 1 / arterial_stiffness

    # Arterial volume from pressure
    arterial_pressure = MAP_internal  # at QSS, pressure = delayed pressure
    arterial_volume = (arterial_pressure - p[23]) * arterial_compliance + p[32]

    # ---- Peripheral and venous volumes/pressures ----
    BB_venous_effect = 1 + p[365] * BB_signal
    venous_compliance = p[40] * p[267] * p[300] * BB_venous_effect

    # Peripheral resistance
    tissue_autoreg_sig1 = -p[107] * (p[105] * (y[27] - p[20] * p[394]) + p[106] * y[42])
    tissue_autoreg = 1 / (0.2 * tissue_autoreg_sig1 + 1)**4
    venous_autoreg_int = 1 - p[109]
    venous_autoreg = venous_autoreg_int + p[109] / (
        1 + np.exp(((tissue_autoreg_sig1 - p[108]) - 1) / p[110]))
    V_ven0_adjusted = p[34] * venous_autoreg

    periph_resist_mult = (p[102] * p[301] * p[302] *
                          (1 - (1 - p[363]) * BB_signal) * tissue_autoreg)
    periph_resist_mult_adj = 1 + p[99] * (periph_resist_mult - 1)
    peripheral_resistance = p[24] * p[100] * periph_resist_mult_adj

    # At QSS, systemic blood flow = CO (conservation)
    # arterial_pressure - peripheral_pressure = CO * peripheral_resistance (roughly)
    # peripheral_pressure - venous_pressure = CO * R_ven0

    # Peripheral volume from pressure balance
    peripheral_pressure = arterial_pressure - CO_m3_per_hr * peripheral_resistance
    peripheral_volume = (peripheral_pressure - p[22]) * p[39] + p[33]

    # Venous volume from blood volume conservation
    blood_volume_m3 = y[33] / 1000  # y[33] is in L, volumes are in m^3

    # Estimate mean LV and RV volumes as roughly initial values
    LV_EDV = y[15]  # use current EDV estimate
    LV_ESV = max(p[48], LV_EDV * 0.4)  # rough ESV
    mean_LV_vol = (LV_EDV + LV_ESV) / 2

    RV_EDV = y[4] if y[4] > 0 else p[123] / p[5]
    mean_RV_vol = RV_EDV

    # Pulmonary volumes from current state (they're relatively stable)
    pulm_art_vol = y[5]
    pulm_ven_vol = y[6]

    if p[416] == 1:  # heart_renal_link enabled
        venous_volume = (blood_volume_m3 - mean_LV_vol - arterial_volume -
                         peripheral_volume - mean_RV_vol - pulm_art_vol - pulm_ven_vol)
        venous_volume = max(venous_volume, p[34] * 0.5)
    else:
        venous_volume = y[0]

    venous_pressure = p[22] + (venous_volume - p[34]) / venous_compliance

    # ---- LV mechanics at end-diastole ----
    LV_fiber_stretch_ed = ((LV_EDV + LV_wall_volume / 3) /
                           (LV_cavity_volume + LV_wall_volume / 3))**0.3333
    LV_sarcomere_length = p[75] * LV_fiber_stretch_ed

    # Passive stress at end-diastole
    stretch_zero_S = p[117] - p[118]
    if LV_fiber_stretch_ed >= stretch_zero_S:
        level_of_hypertrophy = LV_wall_volume / (
            baseline_total_myocyte_volume + total_nonmyocyte_volume)
        hypertrophy_Cf = p[89] * max(0, level_of_hypertrophy - 1)
        C_f = p[65] * (1 + hypertrophy_Cf)
        LV_passive_stress = p[66] * (np.exp(C_f * (
            LV_fiber_stretch_ed - stretch_zero_S)) - 1)
    else:
        LV_passive_stress = 0.0

    LV_radial_stretch = 1 / (LV_fiber_stretch_ed**2)
    if LV_radial_stretch >= 1:
        LV_passive_radial = p[51] * (np.exp(p[50] * (LV_radial_stretch - 1)) - 1)
    else:
        LV_passive_radial = 0.0

    # EDP from passive stress
    if LV_EDV > p[48]:
        rel_vol = 1 + LV_wall_volume / LV_EDV
    else:
        rel_vol = 1 + LV_wall_volume / p[48]
    LV_EDP = (LV_passive_stress - 2 * LV_passive_radial) * np.log(rel_vol) / 3
    LV_EDS = LV_passive_stress

    # Active stress peak (from sarcomere length effect)
    if LV_sarcomere_length > p[73]:
        sl_effect = (LV_sarcomere_length - p[73]) / (p[74] - p[73])
    else:
        sl_effect = 0.0
    beta_contr = 1 - (1 - p[369]) * BB_signal
    peak_active_stress = p[68] * 1.0 * p[67] * sl_effect * 1.0 * beta_contr * p[299]

    # ---- RV mechanics ----
    RV_fiber_stretch = ((y[4] + p[58] / 3) / (p[54] + p[58] / 3))**0.333
    RV_sarcomere_length = p[62] * RV_fiber_stretch

    # ---- Pulmonary pressures at QSS ----
    pulm_art_pressure = (y[5] - p[36]) / p[42] + p[23]
    pulm_ven_pressure = p[22] + (y[6] - p[37]) / (p[41] * BB_venous_effect)

    # ---- Systolic/diastolic pressures at QSS ----
    # At QSS, SBP and DBP are stable. Use stroke volume to estimate pulse pressure.
    SV_m3 = CO_m3_per_hr / heart_rate  # stroke volume
    if arterial_compliance > 0:
        pulse_pressure = SV_m3 / arterial_compliance  # in internal units
    else:
        pulse_pressure = 0.0

    systolic_pressure = arterial_pressure + pulse_pressure * 2 / 3
    diastolic_pressure = arterial_pressure - pulse_pressure / 3

    # Venous pulse pressure (much smaller)
    if venous_compliance > 0:
        ven_pulse = SV_m3 / venous_compliance * 0.1  # much dampened
    else:
        ven_pulse = 0.0
    systolic_venous = venous_pressure + ven_pulse
    diastolic_venous = venous_pressure - ven_pulse

    # ---- Build QSS dict ----
    qss = {
        0: venous_volume,                          # venous_volume
        1: mean_LV_vol,                            # LV_volume (mean)
        2: arterial_volume,                        # arterial_volume
        3: peripheral_volume,                      # peripheral_circulation_volume
        4: mean_RV_vol,                            # RV_volume
        5: pulm_art_vol,                           # pulmonary_arterial_volume
        6: pulm_ven_vol,                           # pulmonary_venous_volume
        7: aortic_flow_m3_hr,                      # aortic_blood_flow_delayed
        8: aortic_flow_m3_hr,                      # pulmonary_blood_flow_delayed
        13: LV_sarcomere_length,                   # LV_sarcomere_length_delayed
        14: RV_sarcomere_length,                   # RV_sarcomere_length_delayed
        15: LV_EDV,                                # LV_EDV
        16: LV_EDP,                                # LV_EDP
        17: LV_EDS,                                # LV_EDS
        18: arterial_pressure,                     # arterial_pressure_delayed
        19: arterial_pressure,                     # arterial_pressure_bigger_delay
        20: systolic_pressure,                     # systolic_pressure
        21: diastolic_pressure,                    # diastolic_pressure
        22: venous_pressure,                       # venous_pressure_delayed
        23: venous_pressure,                       # venous_pressure_bigger_delay
        24: systolic_venous,                       # systolic_venous_pressure
        25: diastolic_venous,                      # diastolic_venous_pressure
    }

    return qss
