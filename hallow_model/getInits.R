getInits <- function() {
  
###Initial conditions - do NOT change order!!!
#Order must match order in model file
#labels are not used by RxODE to match init to compartment
inits <- c(
  venous_volume = as.numeric(theta$venous_volume_0) - 0.0001, ##
  LV_volume = as.numeric(theta$LV_EDV_nom),
  arterial_volume = as.numeric(theta$arterial_volume_0) + 0.0001, ##account for shift of fluid from venous to arterial system when heart is pumping
  peripheral_circulation_volume = theta$peripheral_circulation_volume_0,
  RV_volume = theta$RV_volume_0, 
  pulmonary_arterial_volume = as.numeric(theta$pulmonary_arterial_volume_0) + 0.00005,  ##account for shift of fluid from venous to arterial system when heart is pumping
  pulmonary_venous_volume = as.numeric(theta$pulmonary_venous_volume_0) - 0.00005,
  aortic_blood_flow_delayed = 0, 
  pulmonary_blood_flow_delayed = 0,
  change_in_myocyte_length = 0, change_in_myocyte_diameter = 0,
  LV_active_stress_peak = 50000,
  sim_time = 0, LV_sarcomere_length_delayed = 2.1e-6,
  RV_sarcomere_length_delayed = 2e-6, 
  LV_EDV = as.numeric(theta$LV_EDV_nom), 
  LV_EDP = 1000,
  LV_EDS = 3000,
  arterial_pressure_delayed =as.numeric(theta$nominal_map_setpoint)/0.0075, 
  arterial_pressure_bigger_delay =as.numeric(theta["nominal_map_setpoint"])/0.0075,
  systolic_pressure = (as.numeric(theta$nominal_map_setpoint) + 21)/0.0075, 
  diastolic_pressure = (as.numeric(theta["nominal_map_setpoint"])-10.5)/0.0075,
  venous_pressure_delayed = 917,
  venous_pressure_bigger_delay = 917,
  systolic_venous_pressure = 1050,
  diastolic_venous_pressure = 850,
  AngI=8.164, AngII=5.17,
  AT1_bound_AngII=16.6, AT2_bound_AngII = 5.5, plasma_renin_concentration=17.845,
  blood_volume_L = as.numeric(theta$blood_volume_nom),
  interstitial_fluid_volume=as.numeric(theta$IF_nom),
  sodium_amount= as.numeric(theta$blood_volume_nom)*as.numeric(theta$ref_Na_concentration), 
  IF_sodium_amount= as.numeric(theta$IF_nom)*as.numeric(theta$ref_Na_concentration), 
  stored_sodium = 0,
  tubulo_glomerular_feedback_effect=1,
  normalized_aldosterone_level=1, 
  preafferent_pressure_autoreg_signal=1, 
  glomerular_pressure_autoreg_signal=1, 
  CO_error=0, Na_concentration_error = 0, 
  normalized_vasopressin_concentration_delayed = 1,
  F0_TGF=as.numeric(theta$nom_LoH_Na_outflow), 
  P_bowmans=as.numeric(theta$Pc_pt_s1_mmHg), 
  oncotic_pressure_difference=as.numeric(theta$nom_oncotic_pressure_difference),
  renal_blood_flow_L_min_delayed=as.numeric(theta$nom_renal_blood_flow_L_min),
  SN_macula_densa_Na_flow_delayed = as.numeric(theta$nom_LoH_Na_outflow)/as.numeric(theta$baseline_nephrons),
  rsna_delayed = 1,
  disease_effects_increasing_Kf= 0, disease_effects_decreasing_CD_PN = 0,
  tubular_length_increase=0, tubular_diameter_increase=0, 
  water_out_s1_delayed=3e-8,
  water_out_s2_delayed=1.9e-8,
  water_out_s3_delayed=1.2e-8,
  reabsorbed_urea_cd_delayed =0,#10e-8,
  UGE = 0,
  serum_creatinine = as.numeric(theta$equilibrium_serum_creatinine)*as.numeric(theta$blood_volume_nom), 
  cumNaExcretion = 0, cumWaterExcretion = 0, cumCreatinineExcretion = 0,
  RTg_compensation = 0,
  SGLT2_inhibition_delayed = 1,
  RUGE_delayed = 0,
  CO = as.numeric(theta$CO_nom),
  CO_delayed = as.numeric(theta$CO_nom),
  postglomerular_pressure_delayed = as.numeric(theta$RIHP0),
  postglomerular_pressure_error = 0
  )
  return(inits)
  
}
