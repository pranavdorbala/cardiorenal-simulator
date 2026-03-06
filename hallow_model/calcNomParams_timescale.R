calcNomParams <- function(){

## Parameter File

#Constants and Unit conversions
nL_mL=1e+06
dl_ml=0.01
L_dL=10
L_mL=1000
L_m3=0.001
m3_mL = 1000000
m_mm = 1000 
g_mg=0.001
ng_mg=1e-06
secs_mins=60
min_hr=60
min_sec = 60
hr_day=24
min_day=1440
Pa_mmHg = 0.0075 # (mmHg)
MW_creatinine=113.12
Pi=3.1416
viscosity_length_constant=1.5e-09 
gamma =  1.16667e-5;  #  viscosity of tubular fluid
water_intake_species_scale = 1

#############################
#Cardiac model parameters

#####Cardiovascular System
CO_nom= 5	#L/min
BV = 0.005   # m^3
P_ven0 = 0
P_art0 = 0
R_per0 = 1.27e+008 # (Pa s /m^3) 
R_ven0 = 5e+006 # (Pa s m^-3)
R_art0 = 5e+006 # (Pa s /m^3) 
R_art_pulm = 3e+006 # (Pa s m^-3) 
R_ven_pulm = 6.44e+006 # (Pa s m^-3) 
R_r_atrium = 1000000 # (Pa*s/m3)
R_left_atrium = 1000000 # (Pa*s/m3)
min_flux = 5e-007 # (m3/s)




V_art0 = 0.00045 	# (m^3)
V_per0 = 0.00042 	# (m3)
V_ven0 = 0.0030 	# (m^3) [0.0020] 
V_pulm0 = 0.00015 # (m^3) 
V_pulm_art0 = 0.00004
V_pulm_ven0 = 0.00025

C_art_initial = 1.10e-008 # (m^3 /Pa) 
C_per = 1.0e-008 	# (m^3 /Pa)
C_ven0 = 1.8e-07	# (m^3 /Pa) [2.0e-007]  
					
C_pulm_ven = 1.65e-007 # (m^3/Pa)
C_pulm_art = 2e-8
C_pulm = 1.65e-7
C_lv = 4e-008 #
L_pulm = 60000 
L_art = 60000

LV_V0_baseline = 0.000052	# (m3) [0.000052]
LV_V0_min = 0.00001 		# (m3)
V_w_0 = 0.00012			# (m3) [0.00017]
c_r_LV = 9  
s_r0 = 200 

RV_systolic_time_fraction = 0.3
RV_V0_min = 0.00001 # (m3) 
RV_V0 = 0.000075 # (m3) 
contractility_RV = 1 # 
cf_RV = 8 #   (governs the shape of the passive fiber stress-stretch curve)'
s_f0_RV = 900 # (Pa)(governs the shape of the fiber passive stress-stretch curve)'
V_w_0_RV = 0.0001 # 
c_r_RV = 9 # () 
s_r0_RV = 200 # (Pa) 
sigma_ar_RV = 55000 # (Pa) 
ls_a0_RV = 0.0000015 # (m) 
v0_RV_contraction_velocity_effect_in_RV = 0.00005 # (m/sec)

HR_heart_rate = 70  # (1/min)
cf = 11	# LV cardiac stiffness
s_f0 = 900 # (Pa)(governs the shape of the fiber passive stress-stretch curve)'
sigma_ar = 55000 # (Pa)
contractility =  1# 0.9	## [0.75]
tau_r_LV_twitch_shape = 0.2	# originally = 0.1, changed to improve PV relationship
tau_d_LV_twitch_shape = 0.2
n_r_LV_twitch_shape = 2 # exponent for rise part of LV twitch
n_d_LV_twitch_shape = 4 	# 4 
ls_a0 = 0.0000015 # (m)
ls_ar_sarcomere_length_effect_in_LV = 0.000002 # (m) 
ls_0_passive_LV_sarcomere_length = 0.0000019 # (m)
v0_LV_contraction_velocity_effect_in_LV = 0.000050 # (m/sec) #0.00005
Cv_contraction_velocity_effect_in_LV = 0

##Initial condition for stroke volume
SVnom = CO_nom*1000/HR_heart_rate
LV_EDV_nom = 110/1e6
LV_ESV_nom = (110 - SVnom)/1e6

################################
#### Hypertrophy Parameters ####
Baseline_Myocyte_Number = 3.3e+9
Baseline_Myocyte_Length = 0.000115 # (m)

Baseline_Myocyte_Diameter = 0.000023313467 #(m)
Baseline_Myocyte_Volume = 4.909090909091e-014 # (m3) 

max_myocyte_diameter_increase = 1.5*Baseline_Myocyte_Diameter #1x baseline diameter
max_myocyte_length_increase = Baseline_Myocyte_Length

kD_HYPERTROPHY =  (2*3e-8)/60  # 2e-08 - change in the run file
kL_HYPERTROPHY = (8*2e-9/60)# 1 # 2e-10 - change in the run file
hypertrophy_Cf_slope = 0.39
hypertrophy_contractility_slope = 0#.1

BNP_factor = .0008 #change to manipulate slope of relationship between ln(BNP) and LVEDP

myo_L_scale = 1 # (scaling factor for myocyte lenght contribution to the change in volume)'
myo_D_scale = 0 #(scaling factor for myocyte diameter contribution to the change in volume)'


LV_active_stress_threshhold = 49550 #49200# 54100 #53900	
LV_passive_stress_along_fiber_threshhold = 4500 #4550 #6000#11800


Baseline_Interstitial_Fibrosis = V_w_0*0.02 # (m^3) initial value for a healthy person,  based on Betrami and Anversa 1994, assumes 2% initial interstitial fibrosis'
Baseline_Replacement_Fibrosis = V_w_0*0.02  # (m^3) initial value based on Betrami and Anversa 1994, assumes 2% initial replacement fibrosis)'
Baseline_Interstitial_Tissue = V_w_0*0.22 # (m^3)(data from olivetti and Anversa 1996. This value stays fairly constant across disease states)'


vascular_responsiveness_scale=1
TPR_scale_peripheral_resistance = 1
compliance_scale_arterial_compliance = 1
disease_effect_on_TPR_peripheral_resistance = 1
contractility_scale_LV_active_stress = 1 
c_contr_loss = 1

Kp_CO = 0.1
Ki_CO = 0.001
tissue_autoreg_scale=0.1 #3

venous_autoregulation_signal_Km = 4; #Picked to be large enough so that there is no effect until the autoregulation signal becomes large (uncontrolled)
venous_autoregulation_signal_scale = 0.5;
venous_autoregulation_signal_slope = 0.75;

#parameters governing the saturation of the tissue autoregulation signal
min_tissue_autoreg = 0.4
Vmax_tissue_autoreg = 2
Km_tissue_autoreg = 1.47
hill_tissue_autoreg = 2

#BP effect on stiffness
Stiffness_BP_slope = 0.01 #0.03636364    
C_art_scale = 1
 
stretch_min_LV_passive_stress_along_fiber = 1
stretch_scale_LV_passive_stress_along_fiber = 0


venous_volume_0 = (V_ven0+C_ven0*(((BV)-LV_V0_baseline-V_art0-V_per0-V_ven0-RV_V0-V_pulm0)/(C_art_initial+C_per+C_ven0+C_lv+C_pulm)))
LV_volume_0 = LV_V0_baseline
arterial_volume_0 = (V_art0+C_art_initial*(((BV)-LV_V0_baseline-V_art0-V_per0-V_ven0-RV_V0-V_pulm0)/(C_art_initial+C_per+C_ven0+C_lv+C_pulm)))
peripheral_circulation_volume_0 = (V_per0+C_per*(((BV)-LV_V0_baseline-V_art0-V_per0-V_ven0-RV_V0-V_pulm0)/(C_art_initial+C_per+C_ven0+C_lv+C_pulm)))
RV_volume_0 = ((RV_V0+V_w_0_RV/3)*(((1/cf)*log(((3*(((BV)-LV_V0_baseline-V_art0-V_ven0-RV_V0-V_pulm0)/(C_art_initial+C_ven0+C_lv+C_pulm))/(log(1+(V_w_0/0.000125)))+2*(-(0.193)))/s_f0)+1)+1)^3)-V_w_0_RV/3)
pulmonary_arterial_volume_0 = (V_pulm_art0+C_pulm_art*(((BV)-LV_V0_baseline-V_art0-V_per0-V_ven0-RV_V0-V_pulm_art0-V_pulm_ven0)/(C_art_initial+C_per+C_ven0+C_lv+C_pulm_art+C_pulm_ven)))
pulmonary_venous_volume_0 = (V_pulm_ven0+C_pulm_ven*(((BV)-LV_V0_baseline-V_art0-V_per0-V_ven0-RV_V0-V_pulm_art0-V_pulm_ven0)/(C_art_initial+C_per+C_ven0+C_lv+C_pulm_art+C_pulm_ven)))


######################### Renal Model Parameters

####Systemic parameters
nominal_map_setpoint=85  		#mmHg
IF_nom = 15				#L
blood_volume_nom = 5			#L
Na_intake_rate=100/24/60		#mEq/min  - 100mmol/day or 2300 mg/day
nom_water_intake = 2.1			#L/day
ref_Na_concentration=140 	   	#mEq/L
glucose_concentration = 5.5 	#mmol/L
plasma_albumin_concentration= 35	#mg/ml
plasma_protein_concentration = 7   	#g/dl
plasma_urea = 0
nom_serum_uric_acid_concentration = 5.5 #mg/dl
equilibrium_serum_creatinine=0.92	#mg/dl

###Eventually should link this with P_venous from cardiac model
P_venous=4					#mmHg

####Renal parameters
nom_renal_blood_flow_L_min = 1	#L/min
baseline_nephrons=2e6
nom_Kf=3.9					#nl/min*mmHg
nom_oncotic_pressure_difference= 28 #mmHg
P_renal_vein=4  				#mmHg
nom_oncotic_pressure_peritubular= 28.05	#mmHg		
interstitial_oncotic_pressure = 5	#mmHg

#Renal Vasculature
nom_preafferent_arteriole_resistance= 14 	#mmHg
nom_afferent_diameter=1.65e-5		 	#mmHg
nom_efferent_diameter=1.1e-05 		#mmHg

#Renal Tubules
Dc_pt_nom  = 27e-6			#m
Dc_lh = 17e-6				#m
Dc_dt = 17e-6				#m
Dc_cd = 22e-6				#m

L_pt_s1_nom = 0.005			#m
L_pt_s2_nom = 0.005			#m
L_pt_s3_nom =0.004			#m
L_lh_des = 0.01 				#m
L_lh_asc = 0.01 				#m
L_dct = 0.005				#m	
L_cd = L_lh_des	

tubular_compliance = 0.2	
Pc_pt_s1_mmHg = 20.2#19.4#13.2 #15			#mmHg
Pc_pt_s2_mmHg = 15
Pc_pt_s3_mmHg = 11			#mmHg
Pc_lh_des_mmHg = 8			#mmHg
Pc_lh_asc_mmHg = 7			#mmHg
Pc_dt_mmHg = 6				#mmHg
Pc_cd_mmHg = 5				#mmHg
P_interstitial_mmHg = 5
nominal_pt_na_reabsorption=0.7	#fraction	# can change these 3 to
nominal_loh_na_reabsorption = 0.8	#fraction	# increase BP
nominal_dt_na_reabsorption=0.5	#fraction
LoH_flow_dependence = 0.75

###Renal Glucose reabsorption
nom_glucose_reabs_per_unit_length_s1 = 5.4e-5
nom_glucose_reabs_per_unit_length_s2 = 0
nom_glucose_reabs_per_unit_length_s3 = 2.8e-5
diabetic_adaptation = 1
maximal_RTg_increase = 0.3
T_glucose_RTg = 6000000


glucose_natriuresis_effect_pt = 0
glucose_natriuresis_effect_cd = 0
glucose_diuresis_effect_pt = 0
glucose_diuresis_effect_cd = 0

###Renal urea reabsorption
urea_permeability_PT = 0.5

####Renal albumin seiving
####Proteinuria
nom_glomerular_albumin_sieving_coefficient = 0.00062
SN_albumin_reabsorptive_capacity = 1.4e-6
Emax_seiving = 4
Gamma_seiving = 3
Km_seiving = 25
max_PT_albumin_reabsorption_rate = 0.1
nom_albumin_excretion_rate = 3.5e-9
nom_GP_seiving_damage = 65
c_albumin = 0.0231 #min/nl, from Dean and Lazzara
seiving_inf = 4.25e-4 #from Dean and Lazarra, calculated for seiving coeff =0.00062 when SNGFR = 50 nl/min


####RAAS Pathway parameters
concentration_to_renin_activity_conversion_plasma = 61 
nominal_equilibrium_PRA = 1000 	 	#fmol/ml/hr
nominal_equilibrium_AngI = 7.5 		#fmol/ml
nominal_equilibrium_AngII = 4.75 		#fmol/ml
nominal_renin_half_life = 0.1733		# (hr)
nominal_AngI_half_life = 0.5/60 		#(hr)
nominal_AngII_half_life = 0.66/60 		#(hr)
nominal_AT1_bound_AngII_half_life = 12/60 #hr
nominal_AT2_bound_AngII_half_life = 12/60 #hr
ACE_chymase_fraction = 0.95     		#% of AngI converted by ACE. The rest is converted by chymase
fraction_AT1_bound_AngII = 0.75    		#assume AngII preferentially binds to AT1 vs AT2


########################################################################################
#The following parameters are calculated at equilibrium using the parameters above
########################################################################################

#This pressure is the setpoint that determines the myogenic response of the preafferent vasculature
nom_preafferent_pressure = nominal_map_setpoint - nom_renal_blood_flow_L_min*nom_preafferent_arteriole_resistance;

#This pressure is the setpoint that determines the myogenic response of the afferent vasculature
nom_glomerular_pressure = nom_preafferent_pressure - nom_renal_blood_flow_L_min*(L_m3*viscosity_length_constant/(nom_afferent_diameter^4)/baseline_nephrons);

#This pressure is the setpoint that determines the tubular pressure-natriuresis response 
nom_postglomerular_pressure = nom_preafferent_pressure - nom_renal_blood_flow_L_min*(L_m3*viscosity_length_constant*(1/(nom_afferent_diameter^4)+1/(nom_efferent_diameter^4))/baseline_nephrons);

RIHP0 = 9.32 #nom_postglomerular_pressure 	

# The rate of sodium excretion must equal the rate of sodium intake. Sodium reabsorption rates vary along the tubule, but based on literature
# measurements we have a good, and literature data provides estimates for these rates. However, there is a precise
# rate of sodium reabsorption required to achieve the equilibrium defined by the parameters above.
# Assuming that reabsorption rates are known in all but one segment of the tubule, the exact rate
# of reabsorption of the remaining segment can be calculated. We chose to calculate the CD rate of reabsorpion based on estimates for
# PT, LoH, and DT reabsorption.  
nom_GFR = nom_Kf*(nom_glomerular_pressure - nom_oncotic_pressure_difference - (Pc_pt_s1_mmHg))/nL_mL*baseline_nephrons;
nom_filtered_sodium_load = nom_GFR/L_mL*ref_Na_concentration;

#Reabsorption of Na through SGLT2 and SGLT1 at baseline
nom_filtered_glucose_load = glucose_concentration*nom_GFR/ 1000
nom_glucose_pt_out_s1 = max(0,nom_filtered_glucose_load-nom_glucose_reabs_per_unit_length_s1*L_pt_s1_nom*baseline_nephrons)
nom_glucose_pt_out_s2 = max(0,nom_glucose_pt_out_s1-nom_glucose_reabs_per_unit_length_s2*L_pt_s2_nom*baseline_nephrons)
nom_glucose_pt_out_s3 = max(0,nom_glucose_pt_out_s2-nom_glucose_reabs_per_unit_length_s3*L_pt_s3_nom*baseline_nephrons)
nom_glucose_reabs_per_unit_length_s1

nom_SGTL2_Na_reabs_mmol_s1 = nom_filtered_glucose_load-nom_glucose_pt_out_s1;
nom_SGTL2_Na_reabs_mmol_s2 = nom_glucose_pt_out_s1-nom_glucose_pt_out_s2;
nom_SGTL1_Na_reabs_mmol = 2*(nom_glucose_pt_out_s2-nom_glucose_pt_out_s3);	
nom_total_SGLT_Na_reabs = nom_SGTL2_Na_reabs_mmol_s1+nom_SGTL2_Na_reabs_mmol_s2+nom_SGTL1_Na_reabs_mmol; #mEq/min

nom_SGLT_fractional_na_reabs = nom_total_SGLT_Na_reabs/nom_filtered_sodium_load
#calculate fractional reabsorption of Na through non-SGLT2 mechansims
nominal_pt_na_reabsorption_nonSGLT = nominal_pt_na_reabsorption - nom_SGLT_fractional_na_reabs;

nom_Na_reabs_per_unit_length = -log(1-nominal_pt_na_reabsorption_nonSGLT)/(L_pt_s1_nom+L_pt_s2_nom+L_pt_s3_nom);
nom_Na_pt_s1_reabs = nom_filtered_sodium_load*(1-exp(-nom_Na_reabs_per_unit_length*L_pt_s1_nom));
nom_Na_pt_out_s1 = nom_filtered_sodium_load - nom_Na_pt_s1_reabs - nom_SGTL2_Na_reabs_mmol_s1 ;

nom_Na_pt_s2_reabs = nom_Na_pt_out_s1*(1-exp(-nom_Na_reabs_per_unit_length*L_pt_s2_nom));
nom_Na_pt_out_s2 = nom_Na_pt_out_s1 - nom_Na_pt_s2_reabs - nom_SGTL2_Na_reabs_mmol_s2;

nom_Na_pt_s3_reabs = nom_Na_pt_out_s2*(1-exp(-nom_Na_reabs_per_unit_length*L_pt_s3_nom));
nom_Na_pt_out_s3 = nom_Na_pt_out_s2 - nom_Na_pt_s3_reabs - nom_SGTL1_Na_reabs_mmol;


nom_PT_Na_outflow = nom_Na_pt_out_s3;

nom_Na_in_AscLoH = nom_PT_Na_outflow/baseline_nephrons;
AscLoH_Reab_Rate =(2*nominal_loh_na_reabsorption*nom_Na_in_AscLoH)/L_lh_des; #osmoles reabsorbed per unit length per minute. factor of 2 because osmoles = 2

nom_LoH_Na_outflow = nom_PT_Na_outflow*(1-nominal_loh_na_reabsorption);
nom_DT_Na_outflow = nom_LoH_Na_outflow*(1-nominal_dt_na_reabsorption);
nominal_cd_na_reabsorption = 1-Na_intake_rate/nom_DT_Na_outflow;



#RBF = (MAP - P_venous)/RVR. Given MAP, P_venous, RBF, and preafferent, afferent, and efferent resistances, the remaining peritubular resistance at steady state can be determined
nom_RVR = (nominal_map_setpoint - P_venous)/nom_renal_blood_flow_L_min
nom_peritubular_resistance = nom_RVR - (nom_preafferent_arteriole_resistance + L_m3*viscosity_length_constant*(1/nom_afferent_diameter^4+1/nom_efferent_diameter^4)/baseline_nephrons);

#Calculate the normal amount of sodium reabsorbed per unit surface area of the PT
PT_Na_reab_perUnitSA_0 = (nom_filtered_sodium_load/baseline_nephrons)* nominal_pt_na_reabsorption/(3.14*Dc_pt_nom*(L_pt_s1_nom+L_pt_s2_nom+L_pt_s3_nom))

#Given the values for baseline MAP and CO above, the baseline TPR required to maintain this MAP and CO can be calculated. Since TPR includes renal vascular resistance, the baseline systemic (non-renal) resistance
#can be calculated from this TPR and the values for baseline renal resistances defined above. 
nom_TPR = nominal_map_setpoint/CO_nom

#Calculation of peritubular ultrafiltration coefficient
tubular_reabsorption = nom_GFR/1000 - nom_water_intake*water_intake_species_scale/60/24  #at SS, water excretion equals water intake
#Both RIHP and Kf are unknown, so we can either assume RIHP and calculate Kf, or vice versa. Since RIHP has been measured experimentally,
#it seems better to assume a normal value for RIHP and calculate Kf
nom_peritubular_cap_Kf = - tubular_reabsorption/(nom_postglomerular_pressure - RIHP0 - (nom_oncotic_pressure_peritubular - interstitial_oncotic_pressure))



#Creatinine synthesisrate at equilibrium
creatinine_synthesis_rate  = equilibrium_serum_creatinine * dl_ml * nom_GFR #Units: mg/min



####RAAS Pathway parameters
#Values for half lives and equilibrium concentrations of RAAS peptides available in the literature and 
# defined above to calculate nominal values for other RAAS parameters not available in the literature:
#ACE activity
#Chymase activity
#AT1 receptor binding rate
#AT2 receptor binding rate
#equilibrium AT1_bound_AngII
#These values are then assumed to be fixed unless specified otherwise.
#Calculating these nominal parameter values initially in a separate file is required so that these parameters can then be varied independently in the main model
nominal_equilibrium_PRC = nominal_equilibrium_PRA/concentration_to_renin_activity_conversion_plasma
nominal_AngI_degradation_rate = log(2)/nominal_AngI_half_life #/hr
nominal_AngII_degradation_rate = log(2)/nominal_AngII_half_life #/hr
nominal_AT1_bound_AngII_degradation_rate = log(2)/nominal_AT1_bound_AngII_half_life
nominal_AT2_bound_AngII_degradation_rate = log(2)/nominal_AT2_bound_AngII_half_life
#ACE converts 95% of AngI, chymase converts the rest
nominal_ACE_activity = (ACE_chymase_fraction*(nominal_equilibrium_PRA - nominal_AngI_degradation_rate*nominal_equilibrium_AngI)/nominal_equilibrium_AngI)#Therapy_effect_on_ACE
nominal_chymase_activity = (1-ACE_chymase_fraction)*(nominal_equilibrium_PRA - nominal_AngI_degradation_rate*nominal_equilibrium_AngI)/nominal_equilibrium_AngI
#75% of bound AngII is AT1, the rest is AT2
nominal_AT1_receptor_binding_rate = fraction_AT1_bound_AngII*(nominal_equilibrium_AngI*(nominal_ACE_activity+nominal_chymase_activity)-nominal_AngII_degradation_rate*nominal_equilibrium_AngII)/nominal_equilibrium_AngII
nominal_AT2_receptor_binding_rate = (1-fraction_AT1_bound_AngII)*(nominal_equilibrium_AngI*(nominal_ACE_activity+nominal_chymase_activity)-nominal_AngII_degradation_rate*nominal_equilibrium_AngII)/nominal_equilibrium_AngII
nominal_equilibrium_AT1_bound_AngII = nominal_equilibrium_AngII*nominal_AT1_receptor_binding_rate/nominal_AT1_bound_AngII_degradation_rate
nominal_equilibrium_AT2_bound_AngII = nominal_equilibrium_AngII*nominal_AT2_receptor_binding_rate/nominal_AT2_bound_AngII_degradation_rate

########################################################################################
#The following parameters were determined indirectly from many different literature studies on the response
#various changes in the system (e.g. drug treatments, infusions of peptide, fluid, sodium, etc.....)
########################################################################################

#Effects of AT1-bound AngII on preafferent, afferent, and efferent resistance, and aldosterone secretion
AT1_svr_slope = 0
AT1_preaff_scale = 0.8 #1
AT1_preaff_slope = 16 
AT1_aff_scale=0.8 #0.5
AT1_aff_slope=16
AT1_eff_scale=0.8 #0.1
AT1_eff_slope=16
AT1_PT_scale = 0#.1
AT1_PT_slope = 16
AT1_aldo_slope = 0.02 # 0.05

AT1_aff_EC50 = 1e-9;  #equilibrium AT1-bound AngII in Mol/L
AT1_eff_EC50 = nominal_equilibrium_AT1_bound_AngII*1e-12;  #equilibrium AT1-bound AngII in Mol/L

Emax_AT1_eff = 0
Emax_AT1_aff = 0
AT1_hill = 15

AngII_effect_on_venous_compliance=1;

#Effects of Aldosterone on distal and collecting duct sodium reabsorption
nominal_aldosterone_concentration=85
aldo_DCT_scale=0
aldo_DCT_slope = 0.5
aldo_CD_scale=0.2 #0.35
aldo_CD_slope = 0.5
aldo_renin_slope = -0.4# 2.5#-0.7

#Effects of ANP
normalized_atrial_NP_concentration = 1
nom_ANP = 50 #pg/ml
ANP_aff_scale = 0.2
ANP_aff_slope = 1
ANP_preaff_scale = 0
ANP_preaff_slope = 1
ANP_eff_scale = 0
ANP_eff_slope = 1
anp_CD_scale =-0.1
anp_CD_slope = 2
ANP_effect_on_venous_compliance = 1
LVEDP_ANP_slope = 20   #Maeda 1998 ANP vs. LVEDP

ANP_infused_concentration = 0
ANP_infusion = 0


#Effects of RSNA
renal_sympathetic_nerve_activity = 1
nom_rsna = 1
rsna_preaff_scale = 0.2
rsna_preaff_slope = 0.25
rsna_PT_scale=0
rsna_PT_slope=1
rsna_CD_scale = 0
rsna_CD_slope = 1
rsna_renin_slope=1
rsna_svr_slope = 0
rsna_HR_slope = 0

sna_effect_on_contractility=1;
SNA_effect_on_venous_compliance=1;
B2sna_effect_on_TPR = 1;
A1sna_effect_on_TPR =  1;


#Osmolarity control of vasopressin secretion
Na_controller_gain=0.05
Kp_VP = 2
Ki_VP = 0.005

nom_ADH_urea_permeability = .98
nom_ADH_water_permeability = .98

#Effects of Vasopressin on water intake and reabsorption
nominal_vasopressin_conc=4
water_intake_vasopressin_scale = 0.25#1.5
water_intake_vasopressin_slope = -0.5


#Magnitude and Steepness of tubuloglomerular feedback
S_tubulo_glomerular_feedback=0.7
F_md_scale_tubulo_glomerular_feedback=6
MD_Na_concentration_setpoint = 63.29 #60.8

#Effect of macula densa sodium flow on renin secretion 
md_renin_A = 1
md_renin_tau = 1

#Responsiveness of renal vasculature to regulatory signals
preaff_diameter_range=0.25
afferent_diameter_range=1.2e-05 
efferent_diameter_range=3e-06 
preaff_signal_nonlin_scale=4
afferent_signal_nonlin_scale=4
efferent_signal_nonlin_scale=4

#RAAS pathway (these parameters can be set to different values than used to calculate the equilibrium state above)
AngI_half_life=0.008333 
AngII_half_life=0.011 
AT1_bound_AngII_half_life=0.2 
AT1_PRC_slope=-0.9 
AT1_PRC_yint=0
AT2_bound_AngII_half_life=0.2
concentration_to_renin_activity_conversion_plasma=61
fraction_AT1_bound_AngII=0.75
nominal_ACE_activity=48.9
 nominal_AT1_receptor_binding_rate=12.1
nominal_AT2_receptor_binding_rate=4.0 
nominal_chymase_activity=1.25  
nominal_equilibrium_AT1_bound_AngII=16.63
nominal_equilibrium_PRC=16.4 
renin_half_life=0.1733 

######### Transfer constants for ODEs - determine speed of processes

C_renal_CV_timescale = 60

#These constants govern the calculation of the cardiac cycle and subsequently calculation of CO and MAP
#They should NOT be scaled when the time scale between the renal and cardiac model is scaled. 
C_cycle = 50
C_cycle2 = 100
C_cycle3 = 100
C_co = 0.1
C_co_delay = 0.25
C_map = 0.25

time_step = 1/C_cycle;

#These parameters govern feedbacks between the cardiac and renal model, and ARE scaled when the time scale ebtween the cardiac and renal model is scaled
C_co_error=1
C_vasopressin_delay = 1 #1

#Na and water transfer between blood, IF
Q_water = 1  ##*60
Q_Na = 1  ##*60
Q_Na_store = 0
max_stored_sodium = 500
C_Na_error=1


C_aldo_secretion=100
C_tgf_reset=0
C_md_flow = 0.06 #Time delay between MD sodium flow and renin secretion
C_tgf=1
C_rbf=100
C_serum_creatinine = 1
C_pt_water=1
C_rsna = 100
C_postglomerular_pressure = 1

#################Therapy effects
HCTZ_effect_on_DT_Na_reabs = 1		# Diuretic - Hydrochlorothiazide
HCTZ_effect_on_renin_secretion = 1		# Diuretic - Hydrochlorothiazide

CCB_effect_on_preafferent_resistance = 1	# Calcium Channel Blocker
CCB_effect_on_afferent_resistance = 1	# Calcium Channel Blocker
CCB_effect_on_efferent_resistance = 1	# Calcium Channel Blocker
pct_target_inhibition_MRA = 0			# Mineralocorticoid Receptor Antagonists
pct_target_inhibition_ARB = 0 		# Angiotensin Receptor Blocker
pct_target_inhibition_ACEi = 0		# ACE inhibitor
pct_target_inhibition_DRI = 0					# Direct Renin Inhibitor

ARB_is_on = 0

BB_TPR_effect = 1
BB_cardiac_relax_effect = 0
BB_venous_compliance_effect= 0
BB_preafferent_R_effect = 1			# Beta Blocker
BB_renin_secretion_effect = 1			# Beta Blocker
BB_HR_effect = 1					# Beta Blocker
BB_contractility_effect = 1			# Beta Blocker
BB_is_on = 0            #Turn on Beta Blocker

k_PD = 0.001            #Speed at which drug takes effect (important for heart rate changes w/BB)  



 
#Normalized_aldo_secretion
K_Na_ratio_effect_on_aldo = 1; 

#Renal autoregulation of glomerular pressure and flow
gp_autoreg_scale=0
preaff_autoreg_scale = 0
myogenic_steepness=2
RBF_autoreg_scale = 0
RBF_autoreg_steepness = 1


#Pressure natiuresis effect 
Kp_PN = 1
Kd_PN = 0
Ki_PN = 0#.0005

max_pt_reabs_rate = 0.995
pressure_natriuresis_PT_scale = 0.5#3, 1
pressure_natriuresis_PT_slope = 1

pressure_natriuresis_LoH_scale = 0	#3, 0
pressure_natriuresis_LoH_slope = 1

pressure_natriuresis_DCT_scale = 0	#3, 0
pressure_natriuresis_DCT_slope = 1

max_cd_reabs_rate = 0.995
pressure_natriuresis_CD_scale = 0.5	#3, 1
pressure_natriuresis_CD_slope=1

RBF_CD_scale = 1
RBF_CD_slope = 0.3

#Rate at which the tubular pressure natriuresis mechanism is lost in diabetes (should be zero or negative number)
CD_PN_loss_rate = 0


water_intake_species_scale = 1
CO_species_scale = 1

#Glomerular pressure effect on glomerular hypertrophy
maximal_glom_surface_area_increase = 0.5
T_glomerular_pressure_increases_Kf  = 120000

#PT sodium reabsorption effects on tubular hypertrophy
maximal_tubule_length_increase = 0#.5	
maximal_tubule_diameter_increase = 0#.25	
T_PT_Na_reabs_PT_length = (1e10)
T_PT_Na_reabs_PT_diameter = (1e10)


#Reduce Kf due to glomerulosclerosis
disease_effects_decreasing_Kf = 0

#Disease effects
disease_effect_on_nephrons = 0	

max_s1_Na_reabs = 7.5e-6
max_s2_Na_reabs = 2e-6
max_s3_Na_reabs = 1
max_deltaLoH_reabs=0.75e-6
CD_Na_reabs_threshold = 7e-7  #2.5e-7

#Treatment Parameters
SGLT2_inhibition = 1
SGLT1_inhibition = 1
C_sglt2_delay = 0.1*60
C_ruge = 0.0001*60#.5
Anhe3 = 0
deltaCanp = 0 # To be used for modifying endogenous normalized ANP concentrations
ANP_effect_on_Arterial_Resistance = 0

loop_diuretic_effect = 1

#Aortic Stenosis Parameters
heart_renal_link = 1
aortic_valve_stenosis = 0
R_art_stenosis_factor = 0
stenosis_rate = 0.005

#Mitral Regurgitation Parameters
mitral_regurgitation = 0
mitral_regurgitation_pressure_diff = 1e10  #set really high so no regurgitation
max_mitral_diff = 19000
min_mitral_diff = 16500
k_mitral_diff = 0.05


#Aortic Regurgitation Parameters
aortic_regurgitation = 0
aortic_regurgitation_pressure_diff = 1e10  #set really high so no regurgitation
max_aortic_diff = 7000
min_aortic_diff = 5500
k_aortic_diff = 0.05

t=sort(ls())
param=sapply(t,names)
for (i in 1:length(t)){
  param[i]=get(t[i])
}
param$param=NULL
param = data.frame(param)

return(param)
}




