/*
 * Hallow et al. Cardiorenal Model — Auto-generated C code
 * 
 * Generated from modelfile_commented.R by rxode_to_c.py
 * 70 state variables, 430 parameters
 *
 * Function signature:
 *   void hallow_rhs(double t, const double *y, double *dydt, const double *p)
 *
 * Compile: gcc -O2 -shared -fPIC -o hallow_rhs.so hallow_rhs.c -lm
 */

#include <math.h>
#include <float.h>

/* Safe exp to prevent overflow */
static inline double safe_exp(double x) {
    if (x > 500.0) return exp(500.0);
    if (x < -500.0) return 0.0;
    return exp(x);
}

/* Safe pow: handles negative base with fractional exponent */
static inline double safe_pow(double base, double e) {
    if (base < 0.0 && e != floor(e)) return -pow(-base, e);
    if (base == 0.0 && e <= 0.0) return 0.0;
    return pow(base, e);
}

#define N_STATE 70
#define N_PARAM 430

void hallow_rhs(double t, const double *y, double *dydt, const double *p) {
  /* State variable indices:
   *   y[0] = venous_volume
   *   y[1] = LV_volume
   *   y[2] = arterial_volume
   *   y[3] = peripheral_circulation_volume
   *   y[4] = RV_volume
   *   y[5] = pulmonary_arterial_volume
   *   y[6] = pulmonary_venous_volume
   *   y[7] = aortic_blood_flow_delayed
   *   y[8] = pulmonary_blood_flow_delayed
   *   y[9] = change_in_myocyte_length
   *   y[10] = change_in_myocyte_diameter
   *   y[11] = LV_active_stress_peak
   *   y[12] = sim_time
   *   y[13] = LV_sarcomere_length_delayed
   *   y[14] = RV_sarcomere_length_delayed
   *   y[15] = LV_EDV
   *   y[16] = LV_EDP
   *   y[17] = LV_EDS
   *   y[18] = arterial_pressure_delayed
   *   y[19] = arterial_pressure_bigger_delay
   *   y[20] = systolic_pressure
   *   y[21] = diastolic_pressure
   *   y[22] = venous_pressure_delayed
   *   y[23] = venous_pressure_bigger_delay
   *   y[24] = systolic_venous_pressure
   *   y[25] = diastolic_venous_pressure
   *   y[26] = CO
   *   y[27] = CO_delayed
   *   y[28] = AngI
   *   y[29] = AngII
   *   y[30] = AT1_bound_AngII
   *   y[31] = AT2_bound_AngII
   *   y[32] = plasma_renin_concentration
   *   y[33] = blood_volume_L
   *   y[34] = interstitial_fluid_volume
   *   y[35] = sodium_amount
   *   y[36] = IF_sodium_amount
   *   y[37] = stored_sodium
   *   y[38] = tubulo_glomerular_feedback_effect
   *   y[39] = normalized_aldosterone_level
   *   y[40] = preafferent_pressure_autoreg_signal
   *   y[41] = glomerular_pressure_autoreg_signal
   *   y[42] = CO_error
   *   y[43] = Na_concentration_error
   *   y[44] = normalized_vasopressin_concentration_delayed
   *   y[45] = F0_TGF
   *   y[46] = P_bowmans
   *   y[47] = oncotic_pressure_difference
   *   y[48] = renal_blood_flow_L_min_delayed
   *   y[49] = SN_macula_densa_Na_flow_delayed
   *   y[50] = rsna_delayed
   *   y[51] = disease_effects_increasing_Kf
   *   y[52] = disease_effects_decreasing_CD_PN
   *   y[53] = tubular_length_increase
   *   y[54] = tubular_diameter_increase
   *   y[55] = water_out_s1_delayed
   *   y[56] = water_out_s2_delayed
   *   y[57] = water_out_s3_delayed
   *   y[58] = reabsorbed_urea_cd_delayed
   *   y[59] = UGE
   *   y[60] = serum_creatinine
   *   y[61] = cumNaExcretion
   *   y[62] = cumWaterExcretion
   *   y[63] = cumCreatinineExcretion
   *   y[64] = RTg_compensation
   *   y[65] = SGLT2_inhibition_delayed
   *   y[66] = RUGE_delayed
   *   y[67] = postglomerular_pressure_delayed
   *   y[68] = postglomerular_pressure_error
   *   y[69] = mitral_valve_leak
   */

  /* Parameter indices (see param_map.json for full list):
   *   p[0] = nL_mL
   *   p[1] = dl_ml
   *   p[2] = L_dL
   *   p[3] = L_mL
   *   p[4] = L_m3
   *   p[5] = m3_mL
   *   p[6] = m_mm
   *   p[7] = g_mg
   *   p[8] = ng_mg
   *   p[9] = secs_mins
   *   p[10] = min_hr
   *   p[11] = min_sec
   *   p[12] = hr_day
   *   p[13] = min_day
   *   p[14] = Pa_mmHg
   *   p[15] = MW_creatinine
   *   p[16] = Pi
   *   p[17] = viscosity_length_constant
   *   p[18] = gamma
   *   p[19] = water_intake_species_scale
   *   ... (430 total)
   */

  /* Variables assigned in if/else blocks or shadowing parameters — pre-declared */
  double ANP = 0.0;
  double AscLoH_Reab_Rate = p[230];  /* init from param, overwritten in model */
  double LV_active_stress_peak_old = 0.0;
  double LV_passive_radial_stress = 0.0;
  double LV_passive_stress_along_fiber = 0.0;
  double LV_peak_stress = 0.0;
  double LV_pressure_diastolic_max = 0.0;
  double LV_stress_diastolic_max = 0.0;
  double LV_twitch_shape = 0.0;
  double LV_volume_maximum = 0.0;
  double RV_passive_radial_stress = 0.0;
  double RV_passive_stress_along_fiber = 0.0;
  double RV_twitch_shape = 0.0;
  double R_art = 0.0;
  double SNGFR_nL_min = 0.0;
  double kD_hypertrophy = 0.0;
  double kL_hypertrophy = 0.0;
  double mitral_regurgitation_pressure_diff = p[421];  /* init from param, overwritten in model */
  double mitral_valve_flow_rate = 0.0;
  double mitral_valve_leak_rate = 0.0;
  double nom_glomerular_albumin_sieving_coefficient = p[184];  /* init from param, overwritten in model */
  double normalized_ANP = 0.0;
  double rel_volume = 0.0;
  double rel_volume_LV = 0.0;
  double sarcomere_length_effect_in_LV = 0.0;
  double sarcomere_length_effect_in_RV = 0.0;
  double sin_signal = 0.0;
  double systemic_pressure_maximum = 0.0;
  double systemic_pressure_maximum_1 = 0.0;
  double systemic_pressure_minimum = 0.0;
  double systemic_pressure_minimum_1 = 0.0;
  double systemic_venous_pressure_maximum = 0.0;
  double systemic_venous_pressure_maximum_1 = 0.0;
  double systemic_venous_pressure_minimum = 0.0;
  double systemic_venous_pressure_minimum_1 = 0.0;
  double tubular_reabsorption = p[238];  /* init from param, overwritten in model */
  double venous_volume_target = 0.0;


  //########## Drug Effects #############
  double ARB_signal = p[362]*(1.0-exp(-p[371]*y[12]));
  double BB_signal = p[370]*(1.0-exp(-p[371]*y[12]));
  double BB_venous_effect = (1.0+p[365]*BB_signal);
  double beta_blocker_effect_on_contractility = 1.0-(1.0-p[369] )*BB_signal;
  double beta_blocker_effect_on_heart_rate = (1.0-(1.0-p[368])*BB_signal);

  //########## Disease Effects #############

  //## Aortic stenosis - approximate by increasing afterload resistance
  if (p[417] == 1.0) {
  R_art = p[26]*(1.0+p[418]*(1.0-exp(-p[419]*y[12])));  /* make sure to set init sim_time = 0 */

  //R_art = R_art0*(1+R_art_stenosis_factor*(sim_time/(365*24)));

  } else {
  R_art = p[26];
  }

  //##Allow mitral regurgitation if pressure differential across mitral valve exceeds threshold (set threshold really high to eliminate any regurgitation)
  if (p[420] == 1.0) {
  //decrease threshold for leak gradually over time
  mitral_regurgitation_pressure_diff = p[423] + p[422]*(exp(-p[424]*y[12]));
  } else {
  mitral_regurgitation_pressure_diff = 1e10;  /* set to very large value */
  }


  //########## Heart Rate #############

  double rsna_HR_intercept = 1.0-p[298];

  double rsna_effect_on_HR = p[298] * p[288] + rsna_HR_intercept;

  double heart_rate = p[64] *rsna_effect_on_HR*beta_blocker_effect_on_heart_rate;

  double beat_duration = p[11] / heart_rate ;
  double beat_time = y[12]/beat_duration - floor(y[12]/beat_duration);
  double periods = floor(y[12]/beat_duration);

  //########## Mean Pressures  #############

  double mean_arterial_pressure_MAP = (y[20]/3.0+y[21]*2.0/3.0)*p[14];
  double mean_venous_pressure = (y[24]/3.0+y[25]*2.0/3.0);

  //########## Vascular Autoregulation #############

  //Local tissue autoregulation
  double tissue_autoregulation_sig1 = -p[107]*(p[105]*(y[27] - p[20]*p[394])+p[106]*y[42]);
  double tissue_autoregulation_signal = 1.0/(safe_pow((0.2*tissue_autoregulation_sig1 + 1.0), 4.0));

  //Venous autoregulation
  //designed to only cause venous constriction when CO is decreased and can no longer be maintained by tissue autoregulation
  double venous_autoregulation_signal_int = 1.0-p[109];
  double venous_autoregulation_signal = venous_autoregulation_signal_int + p[109]/(1.0+exp(((tissue_autoregulation_sig1 - p[108])-1.0)/p[110]));
  double V_ven0_adjusted = p[34]*venous_autoregulation_signal;

  //########## Volume unit conversions ###########
  double LV_volume_mL = y[1] * p[5];
  double arterial_volume_mL = y[2] * p[5];
  double peripheral_volume_mL = y[3] * p[5];
  double RV_volume_mL = y[4] * p[5];
  double pulmonary_arterial_volume_mL = y[5] * p[5];
  double venous_volume_mL = y[0] * p[5];
  double total_blood_volume_mL = LV_volume_mL + arterial_volume_mL + peripheral_volume_mL + RV_volume_mL + pulmonary_arterial_volume_mL + y[6] * p[5] + venous_volume_mL;
  double blood_volume = y[33]/1000.0;


  //################################################################################################
  //#################### Cardiac Sub-Model Non-ODE Equations  Part 1 ######################################


  //########## Cardiac tissue composition ###########

  double baseline_total_myocyte_volume = p[49] -                                 p[96] -                                 p[97] -                                 p[98];  /* # baseline myocyte volume determined by V_w_0 */

  double baseline_single_myocyte_volume = baseline_total_myocyte_volume/p[81];

  //baseline myocyte diameter determined by V_w_0 & Baseline_Myocyte_Length, NOT by Baseline_Myocyte_Diameter
  double baseline_myocyte_diameter = 2.0*sqrt(baseline_single_myocyte_volume/(p[16]*p[82]));	

  double myocyte_length = p[82] +                 y[9];		/* # change_in_myocyte_length depends on passive stress levels */

  double myocyte_diameter = baseline_myocyte_diameter +                 y[10];	/* # change_in_myocyte_diameter depends on active stress levels - HYPERTROPHY */

  double single_myocyte_volume = myocyte_length * p[16] * (safe_pow(myocyte_diameter, 2.0)) / 4.0;		/* # myocyte volume calculated as a cylinder */

  double number_of_live_myocytes = p[81];

  double total_myocyte_volume = single_myocyte_volume * number_of_live_myocytes;

  double total_nonmyocyte_volume = p[96] +                           p[98] +                           p[97];

  double LV_wall_volume = total_myocyte_volume +                   total_nonmyocyte_volume;

  //# a measure of how much LV wall wolume has grown
  double level_of_hypertrophy = LV_wall_volume /                       (baseline_total_myocyte_volume + total_nonmyocyte_volume);		

  double pct_change_in_myocyte_diameter = 100.0 * (y[10] / baseline_myocyte_diameter);
  double pct_change_in_myocyte_length = 100.0 * (y[9] / p[82] ) ;


  //########## Cardiac Mechanics ###########
  //# Muscle fiber stress and strain are approximately homogeneously distributed, so that they may be approximated by single values.
  //# Microscopic constitutive laws for fiber stress and radial stress are used to model active and passive fiber stress.

  double LV_cavity_volume = p[47] *                   (safe_pow((1.0 + p[92] * y[9] / p[82]), 3.0)) *                    (safe_pow((1.0 - p[93] * y[10] / baseline_myocyte_diameter), 2.0));

  double LV_fiber_stretch = safe_pow(((y[1] + (LV_wall_volume/3.0)) / (LV_cavity_volume + (LV_wall_volume / 3.0))), 0.3333333);

  double outward_growth = LV_cavity_volume / p[47];		/* # a measure of how much the LV chamber has grown in volume */

  double LV_sarcomere_length = p[75] * LV_fiber_stretch;

  double LV_sarcomere_contraction_velocity = (LV_sarcomere_length - y[13]) / p[336];

  double contraction_velocity_effect_in_LV = (1.0 - LV_sarcomere_contraction_velocity / p[76]) /                                     (1.0 + p[77] * LV_sarcomere_contraction_velocity / p[76]);

  if (LV_sarcomere_length > p[73]) {
  sarcomere_length_effect_in_LV = ((LV_sarcomere_length - p[73]) /                                     (p[74] - p[73]));
  } else {
  sarcomere_length_effect_in_LV = 0.0;
  }

  double chamber_radius = (safe_pow((LV_cavity_volume * 3.0 / 4.0 / p[16]), 0.3333333)) * p[6] ;			/* # approximating the LV as a spherical shell */

  double chamber_diameter = 2.0 * chamber_radius;

  double outer_radius = (safe_pow(((LV_cavity_volume + LV_wall_volume) * 3.0 / 4.0 / p[16]), 0.3333333)) * p[6] ;

  double h_wall = outer_radius - chamber_radius;			/* # wall thickness */

  double h_over_r = h_wall / chamber_radius;				/* # a measure of the LV chamber growth; it's the ratio of wall thickness to chamber radius */

  double EDV_chamber_radius = (safe_pow((y[15]* 3.0 / 4.0 / p[16]), 0.33333333)) * p[6];

  double EDV_chamber_diameter = 2.0*EDV_chamber_radius;

  double EDV_outer_radius = (safe_pow(((y[15] + LV_wall_volume) * 3.0 / 4.0 / p[16]), 0.3333333)) * p[6] ;

  double EDV_h_wall = EDV_outer_radius - EDV_chamber_radius;			/* # wall thickness */

  double EDV_h_over_r = EDV_h_wall / EDV_chamber_radius;	

  double LV_mass = 1000000.0*LV_wall_volume*1.05 ;			/* #wall volume*[(cm->m)^3]*density */

  double LVID = safe_pow(((6.0*y[15])/3.14159), (1.0/3.0));

  //########## Cardiac excitation ###########

  //##Generate justified sinusoidal activation signal - Left Ventricle
  double RV_twitch_duration = p[52] * beat_duration;

  double t_d = p[70]*(1.0+p[364]*BB_signal); 

  double t_r = p[69];

  double t_twitch = t_r + t_d;			

  if (beat_time <= t_r) {
  sin_signal = safe_pow((sin(p[16]*beat_time/t_twitch)), p[71]);
  } else {
  sin_signal = safe_pow((sin(p[16]*beat_time/t_twitch)), p[71]);
  }

  LV_twitch_shape = sin_signal;

  if (beat_time < 0.0) {
  LV_twitch_shape = 0.0;
  }

  if (beat_time > t_twitch) {
  LV_twitch_shape = 0.0;
  }

  //##Generate justified sinusoidal activation signal - Right Ventricle
  RV_twitch_shape = safe_pow((sin(p[16] * beat_time / RV_twitch_duration)), 2.0);

  if (beat_time < 0.0) {
  RV_twitch_shape = 0.0;
  }

  if (beat_time > RV_twitch_duration) {
  RV_twitch_shape = 0.0;
  }

  //########## Left Ventricle Stress Generation ###########

  double LV_active_stress = p[68] *                     LV_twitch_shape *                     p[67] *                     sarcomere_length_effect_in_LV *                     contraction_velocity_effect_in_LV *                     (beta_blocker_effect_on_contractility * p[299]);

  //Increase in left ventricle stiffness with increasing hypertrophy
  double hypertrophy_effect_on_Cf = p[89]*fmax(0.0,(level_of_hypertrophy - 1.0));

  double C_f = p[65]*(1.0+hypertrophy_effect_on_Cf);

  //Unpressurized fiber stretch
  double stretch_zero_S = p[117] -                   p[118];

  if (LV_fiber_stretch >= stretch_zero_S) {
  LV_passive_stress_along_fiber = p[66] *                                   (exp(C_f * (LV_fiber_stretch - stretch_zero_S)) - 1.0);
  } else {
  LV_passive_stress_along_fiber = 0.0;
  }

  double LV_radial_stretch = 1.0/ (LV_fiber_stretch * LV_fiber_stretch);

  if (LV_radial_stretch >= 1.0) {
  LV_passive_radial_stress = p[51] *                               (exp(p[50] * (LV_radial_stretch - 1.0)) - 1.0);
  } else {
  LV_passive_radial_stress = 0.0;
  }

  //Total stress is sum of active and passive stresses
  double LV_total_stress = LV_active_stress +                   LV_passive_stress_along_fiber -                   2.0 * LV_passive_radial_stress;


  if (y[1] > p[48]) {
  rel_volume_LV = 1.0+LV_wall_volume/y[1];
  } else {
  rel_volume_LV = 1.0+LV_wall_volume/p[48];
  }

  //LV Pressure depends on stress and relativ wall volume
  double LV_pressure = LV_total_stress * log(rel_volume_LV)/ 3.0;



  //########## Right Ventricle Stress Generation ###########

  double RV_Cavity_Volume = p[54];

  double RV_wall_volume = p[58];

  double RV_fiber_stretch = safe_pow(((y[4] + p[58]/3.0) / (RV_Cavity_Volume + RV_wall_volume/3.0)), (0.333));

  double RV_sarcomere_length = p[62] * RV_fiber_stretch;

  if (RV_sarcomere_length > p[62]) {
  sarcomere_length_effect_in_RV = (RV_sarcomere_length - p[62]) / (0.000002 - p[62]);
  } else {
  sarcomere_length_effect_in_RV = 0.0;
  }

  double RV_sarcomere_contraction_velocity = (RV_sarcomere_length - y[14]) /                                     p[336];

  double contraction_velocity_effect_in_RV = (1.0 - RV_sarcomere_contraction_velocity / p[63]) /                                     (1.0 + 0.0 * RV_sarcomere_contraction_velocity / p[63]);

  double RV_active_stress_multiplier = p[55]*p[299];

  double RV_active_stress = p[55] *                     RV_twitch_shape *                     p[61] *                     sarcomere_length_effect_in_RV *                     contraction_velocity_effect_in_RV*                     p[299];

  double RV_radial_stretch = 1.0/ (RV_fiber_stretch * RV_fiber_stretch);

  if (RV_radial_stretch >= 1.0) {
  RV_passive_radial_stress = p[60] * (exp(p[59] * (RV_radial_stretch - 1.0)) - 1.0);
  } else {
  RV_passive_radial_stress = 0.0;
  }


  if (RV_fiber_stretch >= 1.0) {
  RV_passive_stress_along_fiber = p[57] * (exp(p[56] * (RV_fiber_stretch - 1.0)) - 1.0);
  } else {
  RV_passive_stress_along_fiber = 0.0;
  }

  double RV_total_stress = RV_active_stress +                   RV_passive_stress_along_fiber -                   2.0 * RV_passive_radial_stress;

  if (y[4] > p[53]) {
  rel_volume = (1.0 + RV_wall_volume / y[4]);
  } else {
  rel_volume = (1.0 + RV_wall_volume / p[53]);
  }

  double RV_pressure = RV_total_stress * log(rel_volume) / 3.0;



  //################################################################################################
  //#################### Vascular Sub-Model Non-ODE Equations  ######################################

  //######## Circulatory hemodynamics #########

  double peripheral_pressure = p[22] +                       (y[3] - p[33]) / p[39];

  double venous_compliance = p[40]*                   p[267]*                   (p[300]*BB_venous_effect);  

  double venous_pressure = p[22] +                   (y[0] - p[34]) / venous_compliance;

  double venous_flow = (peripheral_pressure - venous_pressure) / p[25];


  //#allow blood volume link between heart and kidney to be turned on/off
  if (p[416] == 1.0) {
  venous_volume_target = blood_volume - y[1] - y[2] - y[3] - y[4] - y[5] - y[6];
  } else {
  venous_volume_target = y[0];
  }

  double tricuspid_valve_flow_rate = fmax((venous_pressure - RV_pressure) / p[29],p[31]);

  double pulmonary_arterial_pressure = ( y[5] - p[36] )/(p[42]) +                               p[23];

  double pulmonary_venous_pressure = p[22] +                             (y[6]  - p[37] )/(p[41]*BB_venous_effect);

  double pulmonary_arterial_blood_flow = (pulmonary_arterial_pressure  - pulmonary_venous_pressure )/ p[28] ;

  double dP = RV_pressure -       pulmonary_arterial_pressure;

  double Zn = p[45] +       p[336] * p[27];

  double pulmonary_blood_flow = (y[8] * p[45] + dP * p[336]) / Zn;

  if (pulmonary_venous_pressure > LV_pressure) {
  //Forward flow into ventricle
  mitral_valve_flow_rate = fmax  (  (pulmonary_venous_pressure - LV_pressure)/p[30]  , p[31] ) ;
  } else {

  //If pressure difference less than regurgitation limit, no regurgitation
  if (LV_pressure - pulmonary_venous_pressure < mitral_regurgitation_pressure_diff) {
  mitral_valve_flow_rate = p[31];
  } else {

  //Allow negative flow back out of the ventricle
  mitral_valve_flow_rate = (pulmonary_venous_pressure - LV_pressure)/p[30];
  }
  }

  //Calculate regurgitation rate
  if (mitral_valve_flow_rate < 0.0) {
  mitral_valve_leak_rate = mitral_valve_flow_rate;
  } else {
  mitral_valve_leak_rate = 0.0;
  }

  double pulmonary_valve_flow_rate = fmax(pulmonary_blood_flow,p[31]);

  double peripheral_resistance_multiplier = p[102] *                                     p[301]*                                     p[302]*(1.0-(1.0-p[363])*BB_signal)*                                     tissue_autoregulation_signal;

  double peripheral_resistance_multiplier_adjusted = 1.0+p[99]*(peripheral_resistance_multiplier-1.0);

  double peripheral_resistance = p[24]* p[100] *                         peripheral_resistance_multiplier_adjusted;

  //# arterial_pressure computed without taking account compliance
  //# arterial_compliance = compliance_scale_arterial_compliance * (C_art ) ;
  //# arterial_pressure = (arterial_volume - V_art0) / arterial_compliance + P_art0 ;

  //# arterial_pressure computed taking account compliance
  //# Constants taken from Safar, M. E., et al. Stiffness of carotid artery wall material and blood pressure in humans: application to antihypertensive therapy and stroke prevention. Stroke 31.3 (2000): 782-790.
  //# We assume a linear effect as follows.
  //# BP_effect_on_stiffness = (arterial_pressure - 85)*Stiffness_BP_slope;
  //# Solving the following system


  //When changing compliance, this will allow it to change gradually rather than instantaneously, preventing computational issues
  double C_art = (p[38]-p[38]*p[116])/         exp(y[12]/(24.0/2.0))+p[38]*p[116];

  double Stiffness0 = 1.0/C_art;

  double arterial_stiffness = Stiffness0*                     (1.0+ (mean_arterial_pressure_MAP - p[126])*p[115]);

  double arterial_compliance = 1.0/arterial_stiffness;

  double arterial_pressure = (y[2] - p[32]) / arterial_compliance + p[23] ;

  peripheral_pressure = p[22] +                       (y[3] - p[33]) / p[39];

  double systemic_blood_flow = (arterial_pressure - peripheral_pressure) / peripheral_resistance;


  double dP_1 = LV_pressure - arterial_pressure;

  double Zn_1 = p[46] + R_art * p[336];

  double aortic_blood_flow = (y[7] * p[46] + dP_1 * p[336]) / Zn_1;

  double aortic_valve_flow_rate = fmax(aortic_blood_flow,p[31]);


  //############################################
  //## Relating BNP and NTP to End Diastolic Stress###

  double BNP = exp(p[91]*((y[17]+1736.0)/5.094)+3.14);    
  double NTP = exp((log(BNP)+1.4919)/1.0694);

  double LVEDP_ANP_effect = exp(fmax(0.0, y[16]*0.0075-10.0)/p[285]);

  if (p[416]  == 1.0) {
  if (p[287] == 1.0) {
  ANP = p[286];
  }  else {
  ANP = p[275]*LVEDP_ANP_effect; /* pg/ml */
  }
  normalized_ANP = ANP/p[275];

  } else {
  ANP = p[275];
  normalized_ANP = 1.0;
  }



  //################################################################################################
  //#################### Systolic / Diastolic Calculations  ######################################

  //######## Capture End Diastolic Pressure, Volume, and Stress #########
  if (beat_time >= (1.0 - .01/beat_duration) && beat_time < 1.0) {

  LV_pressure_diastolic_max = LV_pressure;

  LV_stress_diastolic_max = LV_passive_stress_along_fiber;

  LV_volume_maximum = y[1];

  } else {

  LV_pressure_diastolic_max = y[16];

  LV_stress_diastolic_max = y[17];

  LV_volume_maximum = y[15];

  }

  double LV_EDP_old = LV_pressure_diastolic_max;

  double LV_EDS_old = LV_stress_diastolic_max;

  double LV_EDV_old = LV_volume_maximum;


  //######## Capture Systolic and Diastolic Pressures #########

  //# Method to find Systolic and Diastolic blood pressure
  //# use current and last two time steps to find local maxima and minima

  //## SBP and DBP
  if (y[18] < y[19]) {
  systemic_pressure_minimum_1 = y[18];
  } else {
  systemic_pressure_minimum_1 = y[21];
  }

  if (y[18] < arterial_pressure) {
  systemic_pressure_minimum = systemic_pressure_minimum_1;
  } else {
  systemic_pressure_minimum = y[21];
  }

  if (y[18] > y[19]) {
  systemic_pressure_maximum_1 = y[18];
  } else {
  systemic_pressure_maximum_1 = y[20];
  }

  if (y[18] > arterial_pressure) {
  systemic_pressure_maximum = systemic_pressure_maximum_1;
  } else {
  systemic_pressure_maximum = y[20];
  }

  double systolic_pressure_old = systemic_pressure_maximum;

  double diastolic_pressure_old = systemic_pressure_minimum;

  //## Venous Systolic and Diastolic Pressures ###

  if (y[22] < y[23]) {
  systemic_venous_pressure_minimum_1 = y[22];
  } else {
  systemic_venous_pressure_minimum_1 = y[25];
  }

  if (y[22] < venous_pressure) {
  systemic_venous_pressure_minimum = systemic_venous_pressure_minimum_1;
  } else {
  systemic_venous_pressure_minimum = y[25];
  }

  if (y[22] > y[23]) {
  systemic_venous_pressure_maximum_1 = y[22];
  } else {
  systemic_venous_pressure_maximum_1 = y[24];
  }

  if (y[22] > venous_pressure) {
  systemic_venous_pressure_maximum = systemic_venous_pressure_maximum_1;
  } else {
  systemic_venous_pressure_maximum = y[24];
  }

  double systolic_venous_pressure_old = systemic_venous_pressure_maximum;

  double diastolic_venous_pressure_old = systemic_venous_pressure_minimum;


  //## Find Peak systolic stress in the LV ###

  if (beat_time >= (t_r*0.8) && beat_time < (t_r*0.85)) {
  LV_peak_stress = LV_active_stress;
  } else {
  LV_peak_stress = y[11];
  }

  //# fixes a numerical computation problem
  if (LV_active_stress > 1.0) {
  LV_active_stress_peak_old = LV_peak_stress;
  } else {
  LV_active_stress_peak_old = y[11];
  }



  //################################################################################################
  //#################### Cardiac Sub-Model Non-ODE Equations  Part 2 ######################################

  //######## LV Hypertrophy #########

  if (y[11] > p[94]) {
  //# myocyte diameter grows if peak active stress exceeds the threshhold
  kD_hypertrophy = (p[87]*p[329]) *                     fmax(0.0, ((p[85]) - y[10])/                             (p[85]));	
  } else {
  //# regression of hypertrophic growth (negative will come from ODE below)
  kD_hypertrophy = (p[87]*p[329]);		
  }

  if (y[17] > p[95]) {
  //# myocyte length grows if passive stre
  kL_hypertrophy = (p[88]*p[329]) * fmax(0.0, ((p[86]) - y[9])/(p[86]));		/* # progression of growth - myocyte length increase */
  } else {              /* #adding a length threshhold to limit volumetric remodeling */
  kL_hypertrophy = 0.0;		/* # myocyte length does not decrease once stretched */
  }



  //################################################################################################
  //#################### Kidney Sub-Model Non-ODE Equations  Part 2 ######################################

  double number_of_functional_glomeruli = p[140];
  double number_of_functional_tubules = p[140]*(1.0-p[402]);



  //######## Renal Vascular Resistance #########

  //## AT1-bound AngII constricts the preafferent, afferent, and efferent arterioles

  double AT1_preaff_int = 1.0 - p[253]/2.0;
  double AT1_effect_on_preaff = AT1_preaff_int + p[253]/(1.0+exp(-(y[30] - p[250])/p[254]));

  double AT1_aff_int = 1.0 - p[255]/2.0;
  double AT1_effect_on_aff = AT1_aff_int + p[255]/(1.0+exp(-(y[30] - p[250])/p[256]));

  double AT1_eff_int = 1.0 - p[257]/2.0;
  double AT1_effect_on_eff = AT1_eff_int + p[257]/(1.0+exp(-(y[30] - p[250])/p[258]));

  //## ANP may dilate preafferent, afferent, and efferent arteriole

  // ANP_preaff_int = 1 + ANP_preaff_scale/2;
  // ANP_effect_on_preaff = ANP_preaff_int - ANP_preaff_scale/(1+exp(-(normalized_atrial_NP_concentration - 1)/ANP_preaff_slope));

  // ANP_aff_int = 1 + ANP_aff_scale/2;
  // ANP_effect_on_aff = ANP_aff_int - ANP_aff_scale/(1+exp(-(normalized_atrial_NP_concentration - 1)/ANP_aff_slope));

  // ANP_eff_int = 1 + ANP_eff_scale/2;
  // ANP_effect_on_eff= ANP_eff_int - ANP_eff_scale/(1+exp(-(normalized_atrial_NP_concentration - 1)/ANP_eff_slope));


  //## RSNA constricts the preafferent vasculature

  double rsna_preaff_int = 1.0 - p[290]/2.0;

  double rsna_effect_on_preaff = rsna_preaff_int + p[290]/(1.0+exp(-(p[288] - p[289])/p[291]));


  //## Preafferent Resistance
  //The resistance of the arcuate, interlobular arterioles, and other vasculature prior the afferent arterioles is represented by a single resistance - the preafferent arteriole resistance
  //The preafferent arterioles respond myogenically to changes in pressure, and also responds to AT1-bound AngII, RSNA, and ANP
  //The dilation/constriction of the arterioles is limited, and thus the total combined effect of all regulators must saturate

  double preaff_arteriole_signal_multiplier = AT1_effect_on_preaff*                                       y[40]*                                       p[355]*                                       rsna_effect_on_preaff*(1.0-(1.0-p[366])*BB_signal);

  double preaff_arteriole_adjusted_signal_multiplier = (1.0/(1.0+exp(p[319]*(1.0-preaff_arteriole_signal_multiplier)))+0.5);

  double preafferent_arteriole_resistance = p[146]*                                     preaff_arteriole_adjusted_signal_multiplier;


  //## Afferent Arteriole Resistance
  //The afferent arteriole responses the tubuloglomerular feedback (calculated later), as well as to AT1-bound AngII and ANP.
  //It may respond myogenically as well. Some studies suggest the upstream portion responds myogenically while the distal portion responds to TGF. Thus, one could consider the
  //myogenically responsive portion as part of the preafferent resistance.
  //The dilation/constriction of the arterioles is limited, and thus the total combined effect of all regulators must saturate

  double nom_afferent_arteriole_resistance = p[4]*p[17]/                                     (safe_pow(p[147], 4.0));

  double afferent_arteriole_signal_multiplier = y[38] *                                         AT1_effect_on_aff *                                         y[41]*                                         p[356];

  double afferent_arteriole_adjusted_signal_multiplier = (1.0/(1.0+exp(p[320]*(1.0-afferent_arteriole_signal_multiplier)))+0.5);

  double afferent_arteriole_resistance = nom_afferent_arteriole_resistance*                                 afferent_arteriole_adjusted_signal_multiplier;

  //## Efferent Arteriole Resistance
  //The efferent arteriole responses to AT1-bound AngII and ANP.
  //The dilation/constriction of the arterioles is limited, and thus the total combined effect of all regulators must saturate

  double nom_efferent_arteriole_resistance = p[4]*p[17]/                                     (safe_pow(p[148], 4.0));

  double efferent_arteriole_signal_multiplier = AT1_effect_on_eff *                                         p[357];

  double efferent_arteriole_adjusted_signal_multiplier = 1.0/(1.0+exp(p[321]*(1.0-efferent_arteriole_signal_multiplier)))+0.5;

  double efferent_arteriole_resistance = nom_efferent_arteriole_resistance*                                 efferent_arteriole_adjusted_signal_multiplier;


  //## Peritubular Resistance
  //Autoregulation of peritubular resistance allows RBF to be autoregulated separately from GFR
  //This is exploratory for now. By default, this effect is turned off by setting RBF_autoreg_scale to zero

  double RBF_autoreg_int = 1.0 - p[376]/2.0;

  double peritubular_autoreg_signal = RBF_autoreg_int +                               p[376]/(1.0+exp((p[139] - y[48])/p[377]));

  double autoregulated_peritubular_resistance = peritubular_autoreg_signal*                                         p[235];

  //## Renal Vascular Resistance
  double renal_vascular_resistance = preafferent_arteriole_resistance +                             (afferent_arteriole_resistance + efferent_arteriole_resistance) / number_of_functional_glomeruli +                             autoregulated_peritubular_resistance;

  //##Renal blood flow
  double renal_blood_flow_L_min = ((mean_arterial_pressure_MAP - (mean_venous_pressure*0.0075-3.16) )/ renal_vascular_resistance); 

  double renal_blood_flow_ml_hr = renal_blood_flow_L_min * 1000.0 * 60.0;


  //##Renal Vasculature Pressures

  double preafferent_pressure = mean_arterial_pressure_MAP -                         renal_blood_flow_L_min*preafferent_arteriole_resistance;

  double glomerular_pressure = mean_arterial_pressure_MAP  -                       renal_blood_flow_L_min * (preafferent_arteriole_resistance + afferent_arteriole_resistance / number_of_functional_glomeruli);

  double postglomerular_pressure = mean_arterial_pressure_MAP  -                           renal_blood_flow_L_min * (preafferent_arteriole_resistance + (afferent_arteriole_resistance+efferent_arteriole_resistance) / number_of_functional_glomeruli);

  //##Autoregulatory signals for preafferent and afferent resistances

  double preaff_autoreg_int = 1.0 - p[374]/2.0;

  double preafferent_pressure_autoreg_function = preaff_autoreg_int+p[374]/(1.0+exp((p[205] - preafferent_pressure)/p[375]));

  double gp_autoreg_int = 1.0 - p[373]/2.0;

  double glomerular_pressure_autoreg_function = gp_autoreg_int+p[373]/(1.0+exp((p[206] - glomerular_pressure)/p[375]));


  //######## Glomerular Filtration #########

  //Glomerular hypertrophy resulting in increased surface area and thus increased Kf is assumed to occur
  //in response to elevated glomerular pressure. A 2 mmHg buffer is built in (i.e. glomerular pressure must be at least 2 mmHg above normal for hypertrophy to begin
  //The increase in Kf saturates and cannot exceed the fractional increase set by maximal_glom_surface_area_increase

  double GP_effect_increasing_Kf = (p[395] - y[51]) *                           fmax(glomerular_pressure/(p[206]+2.0) - 1.0,0.0) /                           (p[396]/p[329]);

  double glomerular_hydrostatic_conductance_Kf = p[141]*(1.0+y[51]);

  //## Glomerular Fitlration Rate

  double net_filtration_pressure = glomerular_pressure -                           y[47] -                           y[46];

  if (net_filtration_pressure <= 0.0) {
  SNGFR_nL_min = 0.001;
  } else {
  SNGFR_nL_min = glomerular_hydrostatic_conductance_Kf *                     net_filtration_pressure;
  }

  //Unit conversion
  double GFR = (SNGFR_nL_min / 1000.0 / 1000000.0 * number_of_functional_tubules);

  double GFR_ml_min = GFR * 1000.0;

  double filtration_fraction = GFR/renal_blood_flow_L_min;

  double serum_creatinine_concentration = y[60]/y[33];

  double creatinine_clearance_rate = GFR_ml_min * p[1] *                             serum_creatinine_concentration; /* Units: mg/min */

  //## Oncotic pressure ###

  double GPdiff = fmax(0.0, glomerular_pressure - (p[191]));

  double GP_effect_on_Seiving = p[186] * safe_pow(GPdiff, p[187]) / (safe_pow(GPdiff, p[187]) + safe_pow(p[188], p[187]));

  //Dean and Lazzara 2006 - Seiving coefficient decreases as GFR increases
  nom_glomerular_albumin_sieving_coefficient = p[193]/(1.0-(1.0-p[193])*exp(-p[192]*SNGFR_nL_min));

  double glomerular_albumin_sieving_coefficient = nom_glomerular_albumin_sieving_coefficient*(1.0 + GP_effect_on_Seiving);

  double SN_albumin_filtration_rate = p[133]*SNGFR_nL_min*1e-6*glomerular_albumin_sieving_coefficient; /* mg/min */

  double SN_albumin_excretion_rate = fmax(0.0, SN_albumin_filtration_rate - p[185])+p[190];

  double albumin_excretion_rate = SN_albumin_excretion_rate*number_of_functional_tubules;


  //Landis Pappenheimer equation used to calculate oncotic pressure at entrance and exit to glomerulus
  //Oncotic pressure is approximated as varying linearly along the glomerulus. Oncotic pressure in the Bowman's space is zero
  //Thus the average pressure difference is the average of the entrance and exit oncotic pressure
  //We do not consider filtration equilibrium

  double Oncotic_pressure_in = 1.629*p[134]+                       0.2935*(safe_pow(p[134], 2.0));

  double SNRBF_nl_min = 1e6*1000.0*renal_blood_flow_L_min/                 number_of_functional_glomeruli;

  double plasma_protein_concentration_out = (SNRBF_nl_min*p[134]-SN_albumin_filtration_rate)/                                     (SNRBF_nl_min-SNGFR_nL_min);

  double Oncotic_pressure_out = 1.629*plasma_protein_concentration_out+                         0.2935*(safe_pow(plasma_protein_concentration_out, 2.0));

  double oncotic_pressure_avg = (Oncotic_pressure_in+Oncotic_pressure_out)/2.0;


  //######## Plasma sodium concentration and vasopressin secretion #########

  //## Plasma sodium concentration

  double Na_concentration = y[35] /                     y[33];
  double IF_Na_concentration = y[36]/                       y[34];

  double sodium_storate_rate = p[341]*((p[342] - y[37])/p[342])*                       (IF_Na_concentration - p[131]);

  //##Control of vasopressin secretion
  //A proportional-integral controller is used to ensure there is no steady state error in sodium concentration
  //Relative gains of the P and I controller must be chosen carefully.
  //In order to permit a steady-state error, the integral controller can be removed. But care should be given then in choosing the proportional gain

  double Na_water_controller = p[303]*                     (p[304]*(Na_concentration - p[131])+p[305]*y[43]);

  double normalized_vasopressin_concentration = 1.0 + Na_water_controller;

  double vasopressin_concentration = p[308] *                             normalized_vasopressin_concentration;

  //Effect of vasopressin on water intake
  double water_intake_vasopressin_int = 1.0-p[309]/2.0;

  double water_intake = p[19]*               (p[130]/60.0/24.0)*               (water_intake_vasopressin_int + p[309]/(1.0+exp((y[44]-1.0)/p[310])));

  double daily_water_intake = (water_intake * 24.0 * 60.0);


  //######## Tubular Flow and Reabsorption #########

  //Length of tubular segments
  double L_pt_s1 = p[153]*(1.0+y[53]);

  double L_pt_s2 = p[154]*(1.0+y[53]);

  double L_pt_s3 = p[155]*(1.0+y[53]);

  double Dc_pt = p[149]*(1.0+y[54]);

  double L_pt = L_pt_s1+L_pt_s2 + L_pt_s3;

  double SN_filtered_Na_load = (SNGFR_nL_min / 1000.0 / 1000000.0)*                       Na_concentration;

  double filtered_Na_load = SN_filtered_Na_load*                     number_of_functional_tubules;

  //## Regulatory effects on reabsorption

  //##Pressure natriuresis effects

  double pressure_natriuresis_signal = fmax(0.001,                                 1.0+p[378]*(postglomerular_pressure - p[208]) +                                 p[380]*y[68] +                                 p[379]*(postglomerular_pressure - y[67]));

  double pressure_natriuresis_PT_int = 1.0 - p[382]/2.0;

  double pressure_natriuresis_PT_effect = fmax(0.001,                                     pressure_natriuresis_PT_int +                                     p[382] /                                     (1.0 + exp(pressure_natriuresis_signal-1.0))); 

  double pressure_natriuresis_LoH_int = 1.0 - p[384]/2.0;

  double pressure_natriuresis_LoH_effect = fmax(0.001,pressure_natriuresis_LoH_int +                                         p[384] /                                         (1.0 + exp((y[67] - p[208]) / p[385]))); 

  double pressure_natriuresis_DCT_magnitude = fmax(0.0,p[386] );

  double pressure_natriuresis_DCT_int = 1.0 - pressure_natriuresis_DCT_magnitude/2.0;

  double pressure_natriuresis_DCT_effect = fmax(0.001,pressure_natriuresis_DCT_int +                                       pressure_natriuresis_DCT_magnitude /                                       (1.0 + exp((y[67] - p[208]) / p[387]))); 

  double pressure_natriuresis_CD_magnitude = fmax(0.0,p[389] *(1.0+y[52]));

  double pressure_natriuresis_CD_int = 1.0 - pressure_natriuresis_CD_magnitude/2.0;

  double pressure_natriuresis_CD_effect = fmax(0.001,pressure_natriuresis_CD_int +                                       pressure_natriuresis_CD_magnitude /                                       (1.0 + exp(pressure_natriuresis_signal-1.0))); 

  double RBF_CD_int = 1.0 - p[391]/2.0;

  double RBF_CD_effect = fmax(0.001, RBF_CD_int +                       p[391]/                       (1.0+exp((renal_blood_flow_L_min - p[139])/p[392])));

  //##AT1-bound AngII effect on PT reabsorption

  double AT1_PT_int = 1.0 - p[259]/2.0;

  double AT1_effect_on_PT = AT1_PT_int + p[259]/(1.0+exp(-(y[30] - p[250])/p[260]));

  //## RSNA effect on PT and CD sodium reabsorption
  double rsna_PT_int = 1.0 - p[292]/2.0;

  double rsna_effect_on_PT = 1.0;

  double rsna_CD_int = 1.0 - p[294]/2.0;

  double rsna_effect_on_CD = rsna_CD_int + p[294]/(1.0+exp((1.0 - p[288])/p[295]));

  //## Aldosterone effect on distal and collecting duct sodium reabsorption

  double aldosterone_concentration = y[39]*                             p[268]; 

  double Aldo_MR_normalised_effect = y[39]*                               (1.0 - p[358]);

  double aldo_DCT_int = 1.0 - p[269]/2.0;

  double aldo_effect_on_DCT = aldo_DCT_int + p[269]/                           (1.0+exp((1.0 - Aldo_MR_normalised_effect)/p[270]));

  double aldo_CD_int = 1.0 - p[271]/2.0;

  double aldo_effect_on_CD = aldo_CD_int + p[271]/                       (1.0+exp((1.0 - Aldo_MR_normalised_effect)/p[272]));

  //##ANP effect on collecting duct sodium reabsorption
  double anp_CD_int = 1.0 - p[282]/2.0;
  double anp_effect_on_CD = anp_CD_int + p[282]/(1.0+exp((1.0 - normalized_ANP)/p[283]));

  //Effect of SGLT2/NHE3 coupling
  double NHE3inhib = p[412]*y[66];

  double pt_multiplier = AT1_effect_on_PT *                 rsna_effect_on_PT *                 pressure_natriuresis_PT_effect*                 (1.0-NHE3inhib);

  double e_pt_sodreab = fmin(1.0,p[220] * pt_multiplier);

  double e_dct_sodreab = fmin(1.0,p[171] *                       aldo_effect_on_DCT*                       pressure_natriuresis_DCT_effect *                       p[353]); 

  double cd_multiplier = aldo_effect_on_CD*                 rsna_effect_on_CD*                 pressure_natriuresis_CD_effect*                 RBF_CD_effect;

  double e_cd_sodreab = fmin(0.9999,p[233]*cd_multiplier*anp_effect_on_CD);

  //######## Proximal Tubule Reabsorption #########
  //##Glucose Filtration and reabsorption in PT
  //Assume glucose reabsorption depends only on availability of SGLT1/2
  //Assume constant amount of reabsorption per unit length through SGLT2 in convoluted PT
  //Assume constant amount of reabsorption per unit length through SGLT1 in straight/recta PT

  //Chosen so that UGE becomes non-zero for plasma_glucose concentration ~8.5 mmol/l
  double glucose_reabs_per_unit_length_s1 = p[173]*                                     y[65]*                                     (1.0+y[64]);

  double glucose_reabs_per_unit_length_s2 = p[174]*                                     y[65]*                                     (1.0+y[64]);

  double glucose_reabs_per_unit_length_s3 = p[175]*                                     (1.0+y[64])*                                     p[409];

  double SN_filtered_glucose_load = p[132]*SNGFR_nL_min / 1000.0 / 1000000.0;  /* mmol/min */

  double glucose_pt_out_s1 = fmax(0.0,SN_filtered_glucose_load-                             glucose_reabs_per_unit_length_s1*L_pt_s1); /* mmol/min */

  double glucose_pt_out_s2 = fmax(0.0,glucose_pt_out_s1-glucose_reabs_per_unit_length_s2*L_pt_s2); /* mmol/min */

  double glucose_pt_out_s3 = fmax(0.0,glucose_pt_out_s2-glucose_reabs_per_unit_length_s3*L_pt_s3); /* mmol/min */

  double RUGE = glucose_pt_out_s3*number_of_functional_tubules*180.0; /* RUGE in mg/min */

  double excess_glucose_increasing_RTg = (p[177] - y[64]) * fmax(RUGE,0.0) /                                 (p[178]/p[329]);

  double osmotic_natriuresis_effect_pt = 1.0-fmin(1.0,RUGE *p[179]);

  double osmotic_natriuresis_effect_cd = 1.0-fmin(1.0,RUGE *p[180]);

  double osmotic_diuresis_effect_pt = 1.0-fmin(1.0,RUGE *p[181]);

  double osmotic_diuresis_effect_cd = 1.0-fmin(1.0,RUGE *p[182]);

  //## PT Sodium filtration and reabsorption
  // Sodium reabsorbed 1:1 with glucose in S1 and S2
  // Sodium reabsorbed 2:1 with glucose in S3
  // Assume for non-SGLT reabsorption, sodium reabsorbed at a constant RATE along the tubule
  // (represents glomerulotubular balance)

  SN_filtered_Na_load = (SNGFR_nL_min / 1000.0 / 1000000.0)*Na_concentration;  	/* mmol/min */

  double SGTL2_Na_reabs_mmol_s1 = SN_filtered_glucose_load-                           glucose_pt_out_s1;   		/* mmol/min */

  double SGTL2_Na_reabs_mmol_s2 = glucose_pt_out_s1-                           glucose_pt_out_s2;			/* mmol/min */

  double SGTL1_Na_reabs_mmol = 2.0*(glucose_pt_out_s2-glucose_pt_out_s3);			/* mmol/min */

  double total_SGLT_Na_reabs = SGTL2_Na_reabs_mmol_s1 +                       SGTL2_Na_reabs_mmol_s2 +                       SGTL1_Na_reabs_mmol; 	/* mmol/min */

  double Na_reabs_per_unit_length = -log(1.0-e_pt_sodreab)/                             (L_pt_s1+L_pt_s2+L_pt_s3); /* non-SGLT2 reabs	#mmol/min */

  double Na_pt_s1_reabs = fmin(p[403],                     SN_filtered_Na_load*                     (1.0-exp(-Na_reabs_per_unit_length*L_pt_s1)));

  double Na_pt_out_s1 = SN_filtered_Na_load -                 Na_pt_s1_reabs -                 SGTL2_Na_reabs_mmol_s1 ;

  double Na_pt_s2_reabs = fmin(p[404],                     Na_pt_out_s1*                     (1.0-exp(-Na_reabs_per_unit_length*L_pt_s2)));

  double Na_pt_out_s2 = Na_pt_out_s1 -               Na_pt_s2_reabs -               SGTL2_Na_reabs_mmol_s2;

  double Na_pt_s3_reabs = fmin(p[405],                   Na_pt_out_s2*                   (1.0-exp(-Na_reabs_per_unit_length*L_pt_s3)));

  double Na_pt_out_s3 = Na_pt_out_s2 -                 Na_pt_s3_reabs -                 SGTL1_Na_reabs_mmol;

  double PT_Na_reabs_fraction = 1.0-Na_pt_out_s3/                         SN_filtered_Na_load;

  //##PT Urea filtration and reabsorption
  double SN_filtered_urea_load = (SNGFR_nL_min / 1000.0 / 1000000.0)*p[135];

  double urea_out_s1 = SN_filtered_urea_load -               p[183]*               (SN_filtered_urea_load/(0.5*((SNGFR_nL_min / 1000.0 / 1000000.0)+y[55]))-p[135])*               y[55]; /* For now, assuming only reabsorbed at the end */

  double urea_out_s2 = urea_out_s1 -               p[183]*               (urea_out_s1/(0.5*(y[55]+y[56]))-p[135])*               y[56]; /* For now, assuming only reabsorbed at the end */

  double urea_out_s3 = urea_out_s2 -               p[183]*               (urea_out_s2/(0.5*(y[56]+y[57]))-p[135])*               y[57]; /* For now, assuming only reabsorbed at the end */

  double urea_reabsorption_fraction = 1.0-urea_out_s3/SN_filtered_urea_load;


  //##PT Water Reabsorption
  double osmoles_out_s1 = 2.0*Na_pt_out_s1 + glucose_pt_out_s1 + urea_out_s1;

  double water_out_s1 = (((SNGFR_nL_min / 1000.0 / 1000000.0)/                   (2.0*SN_filtered_Na_load+SN_filtered_glucose_load+ SN_filtered_urea_load)))*                 osmoles_out_s1;

  double osmoles_out_s2 = 2.0*Na_pt_out_s2 +                  glucose_pt_out_s2 +                  urea_out_s2;

  double water_out_s2 = (water_out_s1/osmoles_out_s1)*                 osmoles_out_s2;

  double osmoles_out_s3 = 2.0*Na_pt_out_s3 +                   glucose_pt_out_s3 +                   urea_out_s3;

  double water_out_s3 = (water_out_s2/osmoles_out_s2)*                 osmoles_out_s3;

  double PT_water_reabs_fraction = 1.0-water_out_s3/                           (SNGFR_nL_min / 1000.0 / 1000000.0);

  //##Concentrations out of PT
  double Na_concentration_out_s1 = Na_pt_out_s1/water_out_s1;

  double Na_concentration_out_s2 = Na_pt_out_s2/water_out_s2;

  double Na_concentration_out_s3 = Na_pt_out_s3/water_out_s3;

  double glucose_concentration_out_s1 = glucose_pt_out_s1/water_out_s1;

  double glucose_concentration_out_s2 = glucose_pt_out_s2/water_out_s2;

  double glucose_concentration_out_s3 = glucose_pt_out_s3/water_out_s3;

  double urea_concentration_out_s1 = urea_out_s1/water_out_s1;

  double urea_concentration_out_s2 = urea_out_s2/water_out_s2;

  double urea_concentration_out_s3 = urea_out_s3/water_out_s3;

  double osmolality_out_s1 = osmoles_out_s1/water_out_s1;

  double osmolality_out_s2 = osmoles_out_s2/water_out_s2;

  double osmolality_out_s3 = osmoles_out_s3/water_out_s3;

  double PT_Na_outflow = Na_pt_out_s3*number_of_functional_tubules;

  //Tubular sodium reabsorption per unit SA as the driver of tubular hypertrophy
  double PT_Na_reab_perUnitSA = SN_filtered_Na_load*e_pt_sodreab/                         (3.14*Dc_pt*(L_pt_s1+L_pt_s2+L_pt_s3));

  double normalized_PT_reabsorption_density = PT_Na_reab_perUnitSA/p[236];
  double PT_Na_reabs_effect_increasing_tubular_length = 0.0;/* (maximal_tubule_length_increase - tubular_length_increase) * max(normalized_PT_reabsorption_density - 1,0) / (T_PT_Na_reabs_PT_length/C_renal_CV_timescale); */
  double PT_Na_reabs_effect_increasing_tubular_diameter = 0.0;/* (maximal_tubule_diameter_increase - tubular_diameter_increase) * max(normalized_PT_reabsorption_density - 1,0) / (T_PT_Na_reabs_PT_diameter/C_renal_CV_timescale); */


  //#################################### Loop of Henle #########################################

  //####Descending Loop of Henle

  double water_in_DescLoH = water_out_s3; /* L/min */

  double Na_in_DescLoH = Na_pt_out_s3;

  double urea_in_DescLoH = urea_out_s3;

  double glucose_in_DescLoH = glucose_pt_out_s3;

  double osmoles_in_DescLoH = osmoles_out_s3; 

  double Na_concentration_in_DescLoH = Na_concentration_out_s3;

  double Urea_concentration_in_DescLoH = urea_concentration_out_s3;

  double glucose_concentration_in_DescLoH = glucose_concentration_out_s3;

  double osmolality_in_DescLoH = osmoles_out_s3/water_out_s3;

  //No solute reabsorption in descending limb
  double Na_out_DescLoH = Na_in_DescLoH;

  double urea_out_DescLoH = urea_in_DescLoH;

  double glucose_out_DescLoH = glucose_in_DescLoH;

  double osmoles_out_DescLoH = osmoles_in_DescLoH;


  //For LoH, baseline osmoles reabsorbed per unit length is calculated from nominal fractional sodium reabsorption (see baseline parameters file)
  //The rate of reabsorption per unit length may be flow-dependent, and may be modulated by tubular pressure-natriuresis
  // If LoH_flow_dependence = 0, then no flow dependence.

  double deltaLoH_NaFlow = fmin(p[406],p[172]*(Na_out_DescLoH-p[229]));

  AscLoH_Reab_Rate = (2.0*p[170]*(p[229]+deltaLoH_NaFlow)*p[415])/p[156]; /* osmoles reabsorbed per unit length per minute. factor of 2 because osmoles = 2 */

  double effective_AscLoH_Reab_Rate = AscLoH_Reab_Rate*pressure_natriuresis_LoH_effect; /* osmoles reabsorbed per unit length per minute. factor of 2 because osmoles = 2*Na */


  //Min function necesssary to ensure that the LoH does not reabsorb more Na than is delivered to it
  double osmolality_out_DescLoH = osmolality_in_DescLoH*exp(fmin(effective_AscLoH_Reab_Rate*p[156],2.0*Na_in_DescLoH)/(water_in_DescLoH*osmolality_in_DescLoH));

  double water_out_DescLoH = water_in_DescLoH*osmolality_in_DescLoH/osmolality_out_DescLoH;

  double Na_concentration_out_DescLoH = Na_out_DescLoH/water_out_DescLoH;

  double glucose_concentration_out_DescLoH = glucose_out_DescLoH/water_out_DescLoH;

  double urea_concentration_out_DescLoH = urea_out_DescLoH/water_out_DescLoH;


  //####Ascending Loop of Henle

  double Na_in_AscLoH = Na_out_DescLoH;

  double urea_in_AscLoH_before_secretion = urea_out_DescLoH;

  double glucose_in_AscLoH = glucose_out_DescLoH;

  double osmoles_in_AscLoH_before_secretion = osmoles_out_DescLoH;

  double water_in_AscLoH = water_out_DescLoH;


  //Urea Secretion --> Assume all urea reabsorbed and secreted only at tip of loop

  double urea_in_AscLoH = urea_in_AscLoH_before_secretion + y[58];

  double urea_concentration_in_AscLoH = urea_in_AscLoH/water_out_DescLoH;

  double osmoles_in_AscLoH = osmoles_in_AscLoH_before_secretion  + y[58]; 

  double osmolality_in_AscLoH = osmoles_in_AscLoH/water_in_AscLoH;

  //Osmolality descreased due to sodium reabsorption along ascending loop
  //min function necessary so that LoH doesn't reabsorb more sodium than is delivered to it
  double osmolality_out_AscLoH = osmolality_in_AscLoH - fmin(p[156]*effective_AscLoH_Reab_Rate, 2.0*Na_in_DescLoH)*(exp(fmin(p[156]*effective_AscLoH_Reab_Rate, 2.0*Na_in_DescLoH)/(water_in_DescLoH*osmolality_in_DescLoH))/water_in_DescLoH);

  double osmoles_reabsorbed_AscLoH = (osmolality_in_AscLoH - osmolality_out_AscLoH)*water_in_AscLoH;

  double Na_reabsorbed_AscLoH = osmoles_reabsorbed_AscLoH/2.0;

  double Na_out_AscLoH = fmax(0.0,Na_in_AscLoH - Na_reabsorbed_AscLoH);



  //Water, glucose, and urea are not reabsorbed along the ascending limb
  double urea_out_AscLoH = urea_in_AscLoH; /* urea secretion accounted for above */

  double glucose_out_AscLoH = glucose_in_AscLoH;

  double water_out_AscLoH = water_in_AscLoH;

  double osmoles_out_AscLoH = osmolality_out_AscLoH*water_out_AscLoH;


  double Na_concentration_out_AscLoH = Na_out_AscLoH/water_out_AscLoH;

  double glucose_concentration_out_AscLoH = glucose_out_AscLoH/water_out_AscLoH;

  double urea_concentration_out_AscLoH = urea_out_AscLoH/water_out_AscLoH;

  double LoH_reabs_fraction = 1.0-Na_out_AscLoH/Na_in_AscLoH;

  double SN_macula_densa_Na_flow = Na_out_AscLoH;

  double MD_Na_concentration = Na_concentration_out_AscLoH;

  double TGF0_tubulo_glomerular_feedback = 1.0 - p[311]/2.0;

  double tubulo_glomerular_feedback_signal = (TGF0_tubulo_glomerular_feedback + p[311] / (1.0 + exp((p[313] - MD_Na_concentration)/ p[312])));



  //############################Distal Convoluted Tubule #######################

  double water_in_DCT = water_out_AscLoH; 

  double Na_in_DCT = Na_out_AscLoH;

  double urea_in_DCT = urea_out_AscLoH;

  double glucose_in_DCT = glucose_out_AscLoH;

  double osmoles_in_DCT = osmoles_out_AscLoH;

  double Na_concentration_in_DCT = Na_concentration_out_AscLoH; 

  double urea_concentration_in_DCT = urea_concentration_out_AscLoH;

  double glucose_concentration_in_DCT = glucose_concentration_out_AscLoH;

  double osmolality_in_DCT = osmolality_out_AscLoH;


  //Assume only sodium reabsorbed along DCT, no water, glucose, or urea reabsorption
  double urea_out_DCT = urea_in_DCT;

  double glucose_out_DCT = glucose_in_DCT;

  double water_out_DCT = water_in_DCT;

  double urea_concentration_out_DCT = urea_out_DCT/water_out_DCT;

  double glucose_concentration_out_DCT = glucose_out_DCT/water_out_DCT;


  //Assume sodium reabsorption at a constant fraction of delivery
  double R_dct = -log(1.0-e_dct_sodreab)/p[158];

  double Na_out_DCT = Na_in_DCT*exp(-R_dct*p[158]);

  double Na_concentration_out_DCT = Na_out_DCT/water_out_DCT;

  double osmolality_out_DCT = 2.0*Na_concentration_out_DCT + glucose_concentration_out_DescLoH + urea_concentration_in_AscLoH;

  double osmoles_out_DCT = osmolality_out_DCT*water_out_DCT;

  double DCT_Na_reabs_fraction = 1.0-Na_out_DCT/Na_in_DCT; 

  //###############################################Collecting Duct###############################

  double water_in_CD = water_out_DCT;

  double Na_in_CD = Na_out_DCT;

  double urea_in_CD = urea_out_DCT;

  double glucose_in_CD = glucose_out_DCT;

  double osmoles_in_CD = osmoles_out_DCT;

  //Use this to turn off osmotic diuresis effect
  //osmoles_in_CD = osmoles_out_DCT - glucose_in_CD;

  double osmolality_in_CD = osmoles_in_CD/water_in_CD;

  double Na_concentration_in_CD = Na_concentration_out_DCT;

  double urea_concentration_in_CD = urea_concentration_out_DCT;

  double glucose_concentration_in_CD = glucose_concentration_out_DCT;

  osmotic_diuresis_effect_cd = 1.0-fmin(1.0,RUGE *p[182]);

  //###Assume sodium reabsorbed, then water follows
  //### Then urea reabsorbed at end
  //### Then additional water reabsorbed following urea reabsorption

  //Assume sodium reabsorbed at fractional rate eta
  double e_cd_sodreab_adj = e_cd_sodreab*osmotic_natriuresis_effect_cd;

  double R_cd = -log(1.0-e_cd_sodreab_adj)/p[159];

  double Na_reabsorbed_CD = fmin(Na_in_CD*(1.0-exp(-R_cd*p[159])),p[407]);

  double Na_out_CD = Na_in_CD-Na_reabsorbed_CD;

  double CD_Na_reabs_fraction = 1.0-Na_out_CD/Na_in_CD; 

  double ADH_water_permeability_old = fmin(0.99999,fmax(0.0,p[307]*normalized_vasopressin_concentration));

  double ADH_water_permeability = normalized_vasopressin_concentration/(0.15+normalized_vasopressin_concentration);



  //Water reabsorption follows gradient but is regulated by ADH

  double osmoles_out_CD = osmoles_in_CD-2.0*(Na_in_CD - Na_out_CD);

  double osmolality_out_CD_before_osmotic_reabsorption = osmoles_out_CD/water_in_CD;

  double water_reabsorbed_CD = ADH_water_permeability*osmotic_diuresis_effect_cd*water_in_CD*(1.0-osmolality_out_CD_before_osmotic_reabsorption/osmolality_out_DescLoH);

  double water_out_CD = water_in_CD-water_reabsorbed_CD;

  double osmolality_out_CD_after_osmotic_reabsorption = osmoles_out_CD/water_out_CD;

  double glucose_concentration_after_urea_reabsorption = glucose_in_CD/water_out_CD;

  double urine_flow_rate = water_out_CD*number_of_functional_tubules;

  double daily_urine_flow = (urine_flow_rate * 60.0 * 24.0);


  double Na_excretion_via_urine = Na_out_CD*number_of_functional_tubules;

  double Na_balance = p[129] - Na_excretion_via_urine;

  double water_balance = daily_water_intake - daily_urine_flow;

  double FENA = Na_excretion_via_urine/filtered_Na_load;

  double PT_fractional_glucose_reabs = (SN_filtered_glucose_load - glucose_pt_out_s3)/SN_filtered_glucose_load;

  double PT_fractional_Na_reabs = (SN_filtered_Na_load - Na_pt_out_s3)/SN_filtered_Na_load;

  double PT_fractional_urea_reabs = ( SN_filtered_urea_load - urea_out_s3)/SN_filtered_urea_load;

  double PT_fractional_water_reabs = ((SNGFR_nL_min / 1000.0 / 1000000.0) - water_out_s3)/(SNGFR_nL_min / 1000.0 / 1000000.0);

  double LoH_fractional_Na_reabs = (Na_in_DescLoH - Na_out_AscLoH)/Na_in_DescLoH;

  double LoH_fractional_urea_reabs = (urea_in_DescLoH-urea_out_AscLoH)/urea_in_DescLoH;

  double LoH_fractional_water_reabs = (water_in_DescLoH - water_out_AscLoH)/water_in_DescLoH;

  double DCT_fractional_Na_reabs = (Na_in_DCT - Na_out_DCT)/Na_in_DCT;

  double CD_fractional_Na_reabs = (Na_in_CD - Na_out_CD)/Na_in_CD;

  double CD_OM_fractional_water_reabs = (water_in_CD - water_out_CD)/water_in_CD;


  //####################Renal Interstitial Hydrostatic pressure

  //#####RIHP can be approximated from Starling's equation for the peritubular capillaries
  //## Flow out of the capillary = Kf_peritubular*(Peritubular pressure - RIHP - oncotic pressure difference)
  //## Assume that any fluid reabsorbed reenters the capillary.
  //## Therefore, RIHP = Peritubular Pressure - (oncotic pressure in peritubular capillary - interstitial oncotic pressure) + tubular reabsorption/KF
  //Peritubular pressure is assumed to equal postglomerular pressure

  //Oncotic pressure at the entrance to the peritubular circulation equals oncotic pressure at the exit of the glomerular
  double Oncotic_pressure_peritubular_in = Oncotic_pressure_out;

  double plasma_protein_concentration_peritubular_out = (SNRBF_nl_min)*p[134]/(SNRBF_nl_min-urine_flow_rate*1e6*1000.0/number_of_functional_glomeruli);

  double Oncotic_pressure_peritubular_out = 1.629*plasma_protein_concentration_peritubular_out+0.2935*(safe_pow(plasma_protein_concentration_peritubular_out, 2.0));

  double oncotic_pressure_peritubular_avg = (Oncotic_pressure_peritubular_in+Oncotic_pressure_peritubular_out)/2.0;

  //The amount of fluid reabsorbed is the difference between the amount filtered and the amount excreted
  tubular_reabsorption = GFR_ml_min/1000.0 - urine_flow_rate;

  //Renal Interstitial Hydrostatic Pressure
  double RIHP = postglomerular_pressure - (oncotic_pressure_peritubular_avg - p[145]) + tubular_reabsorption/p[239];




  //################ Tubular Pressure #####################

  //####See written documentation for derivation of the equations below
  //flow rates expressed in m3/min, rather than L/min

  double mmHg_Nperm2_conv = 133.32;

  double Pc_pt_s1 = p[161]*mmHg_Nperm2_conv;

  double Pc_pt_s2 = p[162]*mmHg_Nperm2_conv;

  double Pc_pt_s3 = p[163]*mmHg_Nperm2_conv;

  double Pc_lh_des = p[164]*mmHg_Nperm2_conv;

  double Pc_lh_asc = p[165]*mmHg_Nperm2_conv;

  double Pc_dt = p[166]*mmHg_Nperm2_conv;

  double Pc_cd = p[167]*mmHg_Nperm2_conv;

  double P_interstitial = 4.9*mmHg_Nperm2_conv;

  double pi = 3.14;

  //##CD
  double B1 = (4.0*p[160]+1.0)*128.0*p[18]/pi;

  double mean_cd_water_flow = (water_in_CD-water_out_CD)/2.0;

  double B2_cd = (safe_pow(Pc_cd, (4.0*p[160])))/(safe_pow(p[152], 4.0));

  double P_in_cd = safe_pow((safe_pow(0.0, (4.0*p[160]+1.0))+B1*B2_cd*(mean_cd_water_flow/1e3)*p[159]), (1.0/(4.0*p[160]+1.0)));

  double P_in_cd_mmHg = (P_in_cd+P_interstitial)/mmHg_Nperm2_conv;




  //## DCT
  double B2_dt = (safe_pow(Pc_dt, (4.0*p[160])))/(safe_pow(p[151], 4.0));

  double P_in_dt = safe_pow((safe_pow(P_in_cd, (4.0*p[160]+1.0))+B1*B2_dt*(water_in_DCT/1e3)*p[158]), (1.0/(4.0*p[160]+1.0)));

  double P_in_dt_mmHg = (P_in_dt+P_interstitial)/mmHg_Nperm2_conv;



  //## Asc LoH
  double B2_lh_asc = (safe_pow(Pc_lh_asc, (4.0*p[160])))/(safe_pow(p[150], 4.0));

  double P_in_lh_asc = safe_pow((safe_pow(P_in_dt, (4.0*p[160]+1.0))+B1*B2_lh_asc*(water_in_AscLoH/1e3)*p[157]), (1.0/(4.0*p[160]+1.0)));

  double P_in_lh_asc_mmHg = (P_in_lh_asc+P_interstitial)/mmHg_Nperm2_conv;

  //## Desc LoH
  double A_lh_des = effective_AscLoH_Reab_Rate/(water_in_DescLoH*osmolality_in_DescLoH);

  double B2_lh_des = (safe_pow(Pc_lh_des, (4.0*p[160])))*(water_in_DescLoH/1e3)/((safe_pow(p[150], 4.0))*A_lh_des);

  double P_in_lh_des = safe_pow((safe_pow(P_in_lh_asc, (4.0*p[160]+1.0))+B1*B2_lh_des*(1.0-exp(-A_lh_des*p[156]))), (1.0/(4.0*p[160]+1.0)));

  double P_in_lh_des_mmHg = (P_in_lh_des+P_interstitial)/mmHg_Nperm2_conv;

  //## PT

  //Treat urea as if reabsorbed linearly along whole length of PT
  double Rurea = (SN_filtered_urea_load - urea_out_s3)/(L_pt_s1+L_pt_s2+L_pt_s3);

  double urea_in_s2 = SN_filtered_urea_load - Rurea*L_pt_s1;

  double urea_in_s3 = SN_filtered_urea_load - Rurea*(L_pt_s1+L_pt_s2);

  double A_na = Na_reabs_per_unit_length; 

  double flow_integral_s3 = 2.0*(Na_pt_out_s2/A_na)*(1.0-exp(-A_na*L_pt_s3)) - (3.0/2.0)*glucose_pt_out_s2*safe_pow(L_pt_s3, 2.0) + urea_in_s3*L_pt_s3 - (1.0/2.0)*Rurea*(safe_pow(L_pt_s3, 2.0));

  double flow_integral_s2 = 2.0*(Na_pt_out_s1/A_na)*(1.0-exp(-A_na*L_pt_s2)) - (1.0/2.0)*glucose_pt_out_s1*safe_pow(L_pt_s2, 2.0) + urea_in_s2*L_pt_s2 - (1.0/2.0)*Rurea*(safe_pow(L_pt_s2, 2.0));

  double flow_integral_s1 = 2.0*(SN_filtered_Na_load/A_na)*(1.0-exp(-A_na*L_pt_s1)) - (1.0/2.0)*SN_filtered_glucose_load*safe_pow(L_pt_s1, 2.0) + SN_filtered_urea_load*L_pt_s1 - (1.0/2.0)*Rurea*(safe_pow(L_pt_s1, 2.0));


  //S3 segment
  double B2_pt_s3 = (safe_pow(Pc_pt_s3, (4.0*p[160])))/(safe_pow(Dc_pt, 4.0));

  double B3_pt_s3 = (water_out_s2/1e3)/osmoles_out_s2;

  double P_in_pt_s3 = safe_pow((safe_pow(P_in_lh_des, (4.0*p[160]+1.0))+B1*B2_pt_s3*B3_pt_s3*flow_integral_s3), (1.0/(4.0*p[160]+1.0)));

  double P_in_pt_s3_mmHg = (P_in_pt_s3+P_interstitial)/mmHg_Nperm2_conv;


  double B2_pt_s2 = (safe_pow(Pc_pt_s3, (4.0*p[160])))/(safe_pow(Dc_pt, 4.0));

  double B3_pt_s2 = (water_out_s1/1e3)/osmoles_out_s1;

  double P_in_pt_s2 = safe_pow((safe_pow(P_in_pt_s3, (4.0*p[160]+1.0))+B1*B2_pt_s2*B3_pt_s2*flow_integral_s2), (1.0/(4.0*p[160]+1.0)));

  double P_in_pt_s2_mmHg = (P_in_pt_s2+P_interstitial)/mmHg_Nperm2_conv;


  double B2_pt_s1 = (safe_pow(Pc_pt_s1, (4.0*p[160])))/(safe_pow(Dc_pt, 4.0));

  double B3_pt_s1 = (SNGFR_nL_min / 1e12)/(2.0*SN_filtered_Na_load+SN_filtered_glucose_load+ SN_filtered_urea_load);

  double P_in_pt_s1 = safe_pow((safe_pow(P_in_pt_s2, (4.0*p[160]+1.0))+B1*B2_pt_s1*B3_pt_s1*flow_integral_s1), (1.0/(4.0*p[160]+1.0)));

  double P_in_pt_s1_mmHg = (P_in_pt_s1+P_interstitial)/mmHg_Nperm2_conv;




  //###################### Aldosterone and Renin Secretion

  //##Aldosterone is secreted in response to AT1-bound AngII and changes in potassium or sodium concentration

  double AT1_aldo_int = 1.0 - p[261]*p[250];

  double AngII_effect_on_aldo = AT1_aldo_int + p[261]*y[30];

  double N_als = (p[372] * AngII_effect_on_aldo );


  //##Renin is secreted in response to decreases in AT1-bound AngII, decreases in MD sodium flow, or increases in RSNA
  //RSNA effect on renin secretion

  double rsna_renin_intercept = 1.0-p[296];

  double rsna_effect_on_renin_secretion = p[296] * p[288] + rsna_renin_intercept;


  //Macula Densa Sodium flow effect on renin secretion
  //This relationship is known to be non-linear, and md_renin_tau can be calibrated based on data on changes in renin as a functoin of sodium intake

  double md_effect_on_renin_secretion = p[314]*exp(-p[315]*(y[49]*p[140] - p[231]));

  //AT1-bound AngII feedback on renin secretion

  double AT1_bound_AngII_effect_on_PRA = (safe_pow(10.0, (p[325] * log10(y[30] / p[250]) + p[326])));

  //Aldo effect on renin secretion

  double aldo_renin_intercept = 1.0-p[273];

  double aldo_effect_on_renin_secretion = aldo_renin_intercept + p[273]*Aldo_MR_normalised_effect;


  //Plasma renin activity

  double plasma_renin_activity = p[194]* y[32]*(1.0-p[361]);

  //Renin secretion

  double renin_secretion_rate = (log(2.0)/p[328])*p[241]*AT1_bound_AngII_effect_on_PRA*md_effect_on_renin_secretion*p[354]*aldo_effect_on_renin_secretion*(rsna_effect_on_renin_secretion*(1.0-p[367]*BB_signal));

  //RAAS degradation rates

  double renin_degradation_rate = log(2.0)/p[328]; 

  double AngI_degradation_rate = log(2.0)/p[322];

  double AngII_degradation_rate = log(2.0)/p[323];

  double AT1_bound_AngII_degradation_rate = log(2.0)/p[324];

  double AT2_bound_AngII_degradation_rate = log(2.0)/p[327];

  //RAAS rate constants
  double ACE_activity = p[246]*(1.0 - p[360]);

  double chymase_activity = p[247];

  double AT1_receptor_binding_rate = p[248]*(1.0-p[359]*ARB_signal);

  double AT2_receptor_binding_rate = p[249];



  dydt[0] = venous_flow  - tricuspid_valve_flow_rate + p[329]*(venous_volume_target - y[0]);  /* d/dt(venous_volume) */

  dydt[1] = mitral_valve_flow_rate - aortic_valve_flow_rate;  /* d/dt(LV_volume) */

  dydt[2] = (aortic_valve_flow_rate) - (systemic_blood_flow);  /* d/dt(arterial_volume) */

  dydt[3] = systemic_blood_flow - venous_flow;  /* d/dt(peripheral_circulation_volume) */

  dydt[4] = (tricuspid_valve_flow_rate) - (pulmonary_valve_flow_rate);  /* d/dt(RV_volume) */

  dydt[5] = pulmonary_valve_flow_rate - pulmonary_arterial_blood_flow;  /* d/dt(pulmonary_arterial_volume) */

  dydt[6] = pulmonary_arterial_blood_flow - mitral_valve_flow_rate;  /* d/dt(pulmonary_venous_volume) */

  dydt[7] = p[331] * (aortic_blood_flow - y[7]);  /* d/dt(aortic_blood_flow_delayed) */

  dydt[8] = p[331] * (pulmonary_blood_flow - y[8]);  /* d/dt(pulmonary_blood_flow_delayed) */

  dydt[9] = kL_hypertrophy * (y[17] / p[95] - 1.0);  /* d/dt(change_in_myocyte_length) */

  dydt[10] = kD_hypertrophy * (y[11] / p[94] - 1.0);  /* d/dt(change_in_myocyte_diameter) */

  dydt[11] = p[332] *(LV_active_stress_peak_old - y[11]);  /* d/dt(LV_active_stress_peak) */

  dydt[12] = 1.0;  /* d/dt(sim_time) */

  dydt[13] = p[330]* (LV_sarcomere_length - y[13]);  /* d/dt(LV_sarcomere_length_delayed) */

  dydt[14] = p[330]* (RV_sarcomere_length - y[14]);  /* d/dt(RV_sarcomere_length_delayed) */

  dydt[15] = p[331] * (LV_EDV_old - y[15]);  /* d/dt(LV_EDV) */

  dydt[16] = p[331] *(LV_EDP_old - y[16]);  /* d/dt(LV_EDP) */

  dydt[17] = p[331] *(LV_EDS_old - y[17]);  /* d/dt(LV_EDS) */

  dydt[18] = p[331] * (arterial_pressure - y[18]);  /* d/dt(arterial_pressure_delayed) */

  dydt[19] = p[331] * (y[18] - y[19]);  /* d/dt(arterial_pressure_bigger_delay) */

  dydt[20] = p[331] * (systolic_pressure_old - y[20]);  /* d/dt(systolic_pressure) */

  dydt[21] = p[331] * (diastolic_pressure_old - y[21]);  /* d/dt(diastolic_pressure) */

  dydt[22] = p[331] * (venous_pressure - y[22]);  /* d/dt(venous_pressure_delayed) */

  dydt[23] = p[331] * (y[22] - y[23]);  /* d/dt(venous_pressure_bigger_delay) */

  dydt[24] = p[331] * (systolic_venous_pressure_old - y[24]);  /* d/dt(systolic_venous_pressure) */

  dydt[25] = p[331] * (diastolic_venous_pressure_old - y[25]);  /* d/dt(diastolic_venous_pressure) */

  dydt[26] = p[333]*(aortic_valve_flow_rate*60.0/p[4] - y[26]);  /* d/dt(CO) */

  dydt[27] = p[334]*(y[26] - y[27]);  /* d/dt(CO_delayed) */

  //RAAS Pathway
  dydt[28] = plasma_renin_activity - (y[28]) * (chymase_activity + ACE_activity) - (y[28]) * AngI_degradation_rate;   /* d/dt(AngI) */

  dydt[29] = y[28] * (chymase_activity + ACE_activity) - y[29] * AngII_degradation_rate - y[29]*AT1_receptor_binding_rate - y[29]* (AT2_receptor_binding_rate);  /* d/dt(AngII) */

  dydt[30] = y[29] * (AT1_receptor_binding_rate) - AT1_bound_AngII_degradation_rate*y[30];  /* d/dt(AT1_bound_AngII) */

  dydt[31] = y[29] * (AT2_receptor_binding_rate) - AT2_bound_AngII_degradation_rate*y[31];  /* d/dt(AT2_bound_AngII) */

  dydt[32] = renin_secretion_rate - y[32] * renin_degradation_rate;  /* d/dt(plasma_renin_concentration) */


  //Change in Interstitial fluid volume over time is determined by the different between water intake and urine outflow
  dydt[33] = p[329] *(water_intake- urine_flow_rate + p[339]*(Na_concentration - IF_Na_concentration));  /* d/dt(blood_volume_L) */

  dydt[34] = p[329] *p[339]*(IF_Na_concentration - Na_concentration);  /* d/dt(interstitial_fluid_volume) */

  //Change in total body sodium over time is determined by the different between sodium intake and excretion
  dydt[35] = p[329] * (p[129] - Na_excretion_via_urine + p[340]*(IF_Na_concentration - Na_concentration));  /* d/dt(sodium_amount) */

  dydt[36] = p[329] *(p[340]*(Na_concentration - IF_Na_concentration) - sodium_storate_rate);  /* d/dt(IF_sodium_amount) */

  dydt[37] = p[329] *sodium_storate_rate;  /* d/dt(stored_sodium) */

  //These equations serve only to delay the input variable by one timestep. This allows the previous value of the input variable to be used in an equation that appears
  //in the code before the input variable was defined
  dydt[38] = p[329]*(tubulo_glomerular_feedback_signal-y[38]);  /* d/dt(tubulo_glomerular_feedback_effect) */

  dydt[39] = p[329]*p[344] * (N_als-y[39]);  /* d/dt(normalized_aldosterone_level) */

  dydt[40] = p[329]*100.0*(preafferent_pressure_autoreg_function - y[40]);  /* d/dt(preafferent_pressure_autoreg_signal) */

  dydt[41] = 0.0;/* C_glomerular_pressure_autoreg_signal*(glomerular_pressure_autoreg_function - glomerular_pressure_autoreg_signal); */  /* d/dt(glomerular_pressure_autoreg_signal) */

  //Error signals for PI controllers of cardiac output and sodium concentration
  dydt[42] = p[329]*p[337]*(y[27]-p[20]);  /* d/dt(CO_error) */

  dydt[43] = p[329]*p[343]*(Na_concentration - p[131]);  /* d/dt(Na_concentration_error) */

  //This equation allows a delay between the secretion of vasopression and its effect on water intake and tubular water reabsorption
  dydt[44] = p[329]*p[338]*(normalized_vasopressin_concentration - y[44]);  /* d/dt(normalized_vasopressin_concentration_delayed) */

  //TGF resetting. If C_tgf_reset = 0, no TGF resetting occurs. If it is greater than zero, the setpoint will change over time and will eventually
  //come to equal the ambient MD sodium flow rate.
  dydt[45] = p[329]* p[345]*(SN_macula_densa_Na_flow*p[140] - y[45]);  /* d/dt(F0_TGF) */

  //As above, these equations allow a variable to be used in equations that appear in the code before the variable was first defined.
  dydt[46] = p[329]*100.0*(P_in_pt_s1_mmHg - y[46]);  /* d/dt(P_bowmans) */

  dydt[47] = p[329]*100.0*(oncotic_pressure_avg - y[47]);  /* d/dt(oncotic_pressure_difference) */

  dydt[48] = p[329]*p[348]*(renal_blood_flow_L_min - y[48]);  /* d/dt(renal_blood_flow_L_min_delayed) */

  dydt[49] = p[329]*p[346]*( SN_macula_densa_Na_flow - y[49]);  /* d/dt(SN_macula_densa_Na_flow_delayed) */

  dydt[50] = p[329]*p[351]*(p[288] - y[50]);  /* d/dt(rsna_delayed) */

  //##Disease effects (turned off by default)
  //Glomerular hypertrophy
  dydt[51] = GP_effect_increasing_Kf;  /* d/dt(disease_effects_increasing_Kf) */

  //Loss of CD pressure natriuresis response over time
  dydt[52] = p[393];  /* d/dt(disease_effects_decreasing_CD_PN) */

  //Tubular hypertrophy
  dydt[53] = PT_Na_reabs_effect_increasing_tubular_length;  /* d/dt(tubular_length_increase) */

  dydt[54] = PT_Na_reabs_effect_increasing_tubular_diameter;  /* d/dt(tubular_diameter_increase) */


  dydt[55] = p[329]* p[350]*(water_out_s1 - y[55]);  /* d/dt(water_out_s1_delayed) */

  dydt[56] = p[329]*p[350]*(water_out_s2 - y[56]);  /* d/dt(water_out_s2_delayed) */

  dydt[57] = p[329]*p[350]*(water_out_s3 - y[57]);  /* d/dt(water_out_s3_delayed) */

  dydt[58] = 0.0;/* C_pt_water*(reabsorbed_urea_cd - reabsorbed_urea_cd_delayed); */  /* d/dt(reabsorbed_urea_cd_delayed) */

  //Urinary glucose excretion
  dydt[59] = safe_pow(p[329], RUGE); /* mg/hr */  /* d/dt(UGE) */

  //Serum Creatinine
  dydt[60] = p[329]*(p[240] - creatinine_clearance_rate);  /* d/dt(serum_creatinine) */

  dydt[61] = p[329]*Na_excretion_via_urine;  /* d/dt(cumNaExcretion) */

  dydt[62] = safe_pow(p[329], urine_flow_rate);  /* d/dt(cumWaterExcretion) */

  dydt[63] = p[329]*creatinine_clearance_rate;  /* d/dt(cumCreatinineExcretion) */

  dydt[64] = p[329]*excess_glucose_increasing_RTg;  /* d/dt(RTg_compensation) */

  dydt[65] = p[329]*p[410]*(p[408] - y[65]);  /* d/dt(SGLT2_inhibition_delayed) */

  dydt[66] = p[329]*p[411]*(RUGE - y[66]);  /* d/dt(RUGE_delayed) */

  dydt[67] = p[329]*p[352]*(postglomerular_pressure - y[67]); /* necessary to prevent large oscillations */  /* d/dt(postglomerular_pressure_delayed) */

  dydt[68] = p[329]*(postglomerular_pressure - p[208]);  /* d/dt(postglomerular_pressure_error) */

  dydt[69] = mitral_valve_leak_rate;  /* d/dt(mitral_valve_leak) */

}

/* Wrapper for scipy/ctypes: void rhs(int *n, double *t, double *y, double *dydt, double *p, int *ip) */
void rhs_desolve(int *n, double *t, double *y, double *dydt, double *rpar, int *ipar) {
    hallow_rhs(*t, y, dydt, rpar);
}
