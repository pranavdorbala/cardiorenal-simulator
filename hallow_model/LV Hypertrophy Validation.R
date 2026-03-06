#Load Libraries
library(RxODE)
library(lattice)

####Load model structure
setwd("C:/models/cardiorenalmodel/test/LV Hypertrophy")
source("modelfile_commented.R")

#Compile model
cvrsim <- RxODE(model = ode)

#Load Helper files
source("calcEndTime.R")
source("calcOutputs.R")

#load basecase parameters
source("calcNomParams_timescale.R")
theta=calcNomParams()

#load basecase parameters
source("getInits.R")
inits=getInits()


########################################## Turn off all feedbacks ########################################
##########################################################################################################

theta["tissue_autoreg_scale"] = 0
theta["pressure_natriuresis_PT_scale"] = 0
theta["pressure_natriuresis_CD_scale"] = 0
theta["pressure_natriuresis_LoH_scale"] = 0
theta["pressure_natriuresis_DCT_scale"] = 0

theta["AT1_preaff_scale"] = 0
theta["AT1_aff_scale"]=0
theta["AT1_eff_scale"]=0
theta["AT1_PT_scale"] = 0
theta["AT1_aldo_slope"] = 0
theta["aldo_DCT_scale"]=0
theta["aldo_CD_scale"]=0
theta["ANP_aff_scale"] = 0
theta["ANP_preaff_scale"] = 0
theta["ANP_eff_scale"] = 0
theta["anp_CD_scale"] =0
theta["rsna_preaff_scale"] = 0
theta["rsna_CD_scale"] = 0
theta["rsna_PT_scale"] = 0
theta["rsna_renin_slope"] = 0
theta["rsna_svr_slope"] = 0
theta["kD_HYPERTROPHY"] = 0
theta["kL_HYPERTROPHY"] = 0
theta["S_tubulo_glomerular_feedback"] = 0
theta["Stiffness_BP_slope"] = 0
theta["Ki_PN"] = 0

########################################## Run to initial steady-state with no feedbacks ####################################
##########################################################################################################



endtime = calcEndTime(as.numeric(theta["HR_heart_rate"]), 24*21)
t = seq(0, endtime , .01)
ev = eventTable()
ev$add.sampling(t)

inits = as.list(inits)
x <- cvrsim$run(theta, ev, inits=inits, atol = 1e-6, rtol = 1e-4)
#Calculate outputs not directly calculated in model file
x = calcOutputs(x)



########################################## Check #1#######################################################
#1) Are min/max or average values of continuous signals being captured properly?
##########################################################################################################

#only plot final values
plottimes = 1:5000
par(mfrow=c(2,2))

plot(x[plottimes ,"time"], x[plottimes ,"LV_volume"]*1e6, type="l", main = "LV EDV", ylab = "mL", xlab="time (hrs)")
points(x[plottimes ,"time"],x[plottimes ,"LV_EDV"]*1e6, col="red", type="l", lwd=3)
points(x[plottimes ,"time"],x[plottimes ,"LV_ESV"]*1e6, col="red", type="l", lwd=3)

plot(x[plottimes ,"time"],x[plottimes ,"LV_pressure"]*0.0075, type="l", main = "LV EDP", ylab="mmHg", xlab="time (hrs)")
points(x[plottimes ,"time"],x[plottimes ,"LV_EDP"]*0.0075, col="red", type="l", lwd=3)

plot(x[plottimes ,"time"],x[plottimes ,"arterial_pressure"]*0.0075, type="l", main = "Mean Arterial Pressure", ylab = "mmHg", xlab="time (hrs)")
points(x[plottimes ,"time"],x[plottimes ,"mean_arterial_pressure_MAP"], col="red", type="l", lwd=3)

plot(x[plottimes ,"time"],x[plottimes ,"aortic_valve_flow_rate"]*60/0.001, type="l", main = "Cardiac Output", ylab = "L/min", xlab="time (hrs)")
points(x[plottimes ,"time"],x[plottimes ,"CO"], col="red", type="l", lwd=3)


########################################## Check #2 #######################################################
#Some variables have parameters that specify their setpoint under normal conditions. Are these variables at
#or very close to their setpoints?
##########################################################################################################


compareWithNoms <- function(varname, nomval, range){
  nomval = as.numeric(theta[nomval])
  vallims= c(nomval*(1-range), nomval*(1+range))
  minmax = c(nomval*(1-5*range), nomval*(1+5*range))
  
  plot(x[,"time"]/24,x[,varname], type="l", main = varname, ylim = minmax, lwd=2)
  polygon(c(0, max(x[,"time"]), max(x[,"time"]),0)/24, c(vallims[1], vallims[1], vallims[2], vallims[2]), col=rgb(1, 0, 0, 0.2), border= NA)
  points(c(0, max(x[,"time"]))/24, c(nomval, nomval), type="l", col="gray", lty = 3, lwd=3)
}

#dev.new()
par(mfrow=c(3,3))
par(mar=c(3,3,3,1))
compareWithNoms("CO_delayed", "CO_nom", 0.01)
compareWithNoms("mean_arterial_pressure_MAP", "nominal_map_setpoint", 0.025)
compareWithNoms("blood_volume_L", "blood_volume_nom", 0.025)
compareWithNoms("interstitial_fluid_volume", "IF_nom", 0.025)
compareWithNoms("renal_blood_flow_L_min", "nom_renal_blood_flow_L_min", 0.025)
compareWithNoms("GFR_ml_min", "nom_GFR", 0.025)
compareWithNoms("Na_concentration","ref_Na_concentration", 0.01)
compareWithNoms("AT1_bound_AngII", "nominal_equilibrium_AT1_bound_AngII", 0.025)
compareWithNoms("aldosterone_concentration", "nominal_aldosterone_concentration", 0.025)


########################################## Check #3 #######################################################
#Some variables do not have setpoints, but have a known physiological range in healthy subjects
#Are these variables within that normal range?
##########################################################################################################


compareWithRange <- function(A, i, range){
  
  vallims= c(A$Min[i], A$Max[i])
  minmax = c(A$Min[i]*(1-range), A$Max[i]*(1+range))
  varname = as.character(A$Variable[i])
  plot(x[,"time"]/24,x[,varname]*A$UnitScale[i], type="l", main = varname, ylim = minmax, lwd=2,ylab=A$Units[i], xlab="Time (days)")
  polygon(c(0, max(x[,"time"]), max(x[,"time"]),0)/24, c(vallims[1], vallims[1], vallims[2], vallims[2]), col=rgb(1, 0, 0, 0.2), border= NA)
}


A = read.table("VariableCheckFile.csv",header=T, sep=",")
Nplots = dim(A)[1]
ncols = 6
nrows = ceiling(Nplots/ncols)

#dev.new()
par(mfrow=c(nrows, ncols))

for (i in 1:Nplots) {
  compareWithRange(A,i,0.125)
  
}


########################################## Check #4 #######################################################
#When feedbacks are turned on, under the baseline steady state conditions, these feedbacks should all be at
#their operating point. They should only move away from zero when the system is perturbed from equilibrium
#Are all feedbacks at their operating point?
##########################################################################################################

compareWithOne <- function(varname, range){
  nomval = 1
  vallims= c(nomval*(1-range), nomval*(1+range))
  minmax = c(nomval*(1-5*range), nomval*(1+5*range))
  
  plot(x[,"time"]/24,x[,varname], type="l", main = "", ylim = minmax, lwd=2, xlab = "time", ylab = "normalized")
  polygon(c(0, max(x[,"time"]), max(x[,"time"]),0)/24, c(vallims[1], vallims[1], vallims[2], vallims[2]), col=rgb(1, 0, 0, 0.2), border= NA)
  points(c(0, max(x[,"time"]))/24, c(nomval, nomval), type="l", col="gray", lty = 3, lwd=3)
}


compareWithZero <- function(varname, range){
  nomval = 0
  vallims= c(0-range, 0+range)
  minmax = c(0-5*range, 0+5*range)
  
  plot(x[,"time"]/24,x[,varname], type="l", main = "", ylim = minmax, lwd=2, xlab = "time", ylab = "microns")
  polygon(c(0, max(x[,"time"]), max(x[,"time"]),0)/24, c(vallims[1], vallims[1], vallims[2], vallims[2]), col=rgb(1, 0, 0, 0.2), border= NA)
  points(c(0, max(x[,"time"]))/24, c(nomval, nomval), type="l", col="gray", lty = 3, lwd=3)
}


####Turn feedbacks back on by reloading parameters. All should now be at 1. 
#Reload parameters
theta=calcNomParams()
theta_orig = theta

endtime = calcEndTime(as.numeric(theta["HR_heart_rate"]), 24*14)
t = seq(0, endtime , .01)
ev = eventTable()
ev$add.sampling(t)


#Save initial conditions for model at steadystate with no feedbacks. Then run with feedbacks turned on.
x = data.frame(x)
inits = as.list(x[dim(x)[1],names(x) %in% names(inits)])

inits["CO_error"] = 0
inits["sim_time"] = 0
inits["postglomerular_pressure_error"] = 0

#Keep remodeling turned off

x <- cvrsim$run(theta, ev, inits=inits)
x = calcOutputs(x)
 x_orig = x
# 
# #save last initial conditions and start there
inits = as.list(x[dim(x)[1],names(x) %in% names(inits)])
inits_orig = inits
theta_orig = theta

dev.new()
#pdf(file = "FigS1.pdf",width = 12, height = 9 )
par(mfrow=c(3,3))
compareWithOne("tissue_autoregulation_signal", 0.05)
title("CO - TPR controller")
compareWithOne("tubulo_glomerular_feedback_effect", 0.05)
title("Tubuloglomerular \n Feedback Signal")
compareWithOne("preafferent_pressure_autoreg_signal", 0.05)
title("Myogenic autoregulation \n signal")
compareWithOne("pressure_natriuresis_PT_effect", 0.05)
title("Pressure Natriuresis Signal")
compareWithOne("normalized_aldosterone_level", 0.05)
title("Vasopressin signal")
compareWithOne("AT1_bound_AngII_effect_on_PRA", 0.05)
title("At1-bound AngII \nfeedback on PRA")
compareWithOne("aldo_effect_on_CD", 0.05)
title("Aldosterone effect \non Na+ reabsorption")

compareWithZero("change_in_myocyte_diameter", 2e-7)
title("Myocyte hypertrophy\nDiameter change")
compareWithZero("change_in_myocyte_length", 2e-7)
title("Myocyte Hypertrophy\nLength Change")

########################################## Check #5 #######################################################
#When the system is perturbed, do cardiac output and sodium balance stabilize reasonably quickly and 
#with sufficient damping?
##########################################################################################################

endtime = calcEndTime(as.numeric(theta["HR_heart_rate"]), 24*28)

t = seq(0, 24, .01)
t = c( t, seq(24.1, endtime, .1))
ev = eventTable()
ev$add.sampling(t)

#Increase cardiac output setpoint from 5 to 5.5 L and test response
theta=calcNomParams()

theta["CO_nom"] = 5.5
x3 <- cvrsim$run(theta, ev, inits=inits)

#dev.new()
par(mfrow=c(1,2))
plot(x3[,"time"]/24,x3[,"CO_delayed"], type="l", main = "Cardiac Output", ylim = c(4.5, 6), lwd=2, ylab = "L/min", xlab = "Day")
polygon(c(0, max(x3[,"time"]), max(x3[,"time"]),0)/24, c(5.4,5.4,5.6,5.6), col=rgb(1, 0, 0, 0.2), border= NA)
points(c(0, max(x3[,"time"]))/24, c(5.5, 5.5), type="l", col="gray", lty = 3, lwd=3)

inits_last = inits

#Increase sodium intake rate and test response
inits = inits_last

theta=calcNomParams()


theta["Na_intake_rate"] = 0.14


t = seq(0, 18*60 , .1)
ev = eventTable()
ev$add.sampling(t)
theta["C_renal_CV_timescale"] = 1

x4_1 <- data.frame(cvrsim$run(theta, ev, inits=inits))
plot(x4_1[,"mean_arterial_pressure_MAP"])
points(xold[,"mean_arterial_pressure_MAP"], col="red")
points(xold2[,"mean_arterial_pressure_MAP"], col = "blue")

inits = as.list(x4_1[dim(x4_1)[1],names(x4_1) %in% names(inits)])
t = seq(0, 24*5, .01)
ev = eventTable()
ev$add.sampling(t)
theta["C_renal_CV_timescale"] = 60
x4_11 <- cvrsim$run(theta, ev, inits=inits)


x4_1[,"time"] = x4_1[,"time"]/60
x4_11[,"time"] = x4_11[,"time"] + max(x4_1[,"time"])
x4 = rbind(x4_1, x4_11)
plot(x4[,"time"],x4[,"mean_arterial_pressure_MAP"])


dev.new()
par(mfrow = c(3,3))
plot(x4[,"time"]/24,x4[,"Na_balance"], type="l", main = "Na_balance", ylim = c(-0.1, 0.1), lwd=2, ylab = "meq/min", xlab = "Day")
polygon(c(0, max(x4[,"time"]), max(x4[,"time"]),0)/24, c(-0.02,-0.02, 0.02, 0.02), col=rgb(1, 0, 0, 0.2), border= NA)
points(c(0, max(x4[,"time"]))/24, c(0, 0), type="l", col="gray", lty = 3, lwd=3)

plot(x4[,"time"]/24,x4[,"mean_arterial_pressure_MAP"], type="l") #, main = "Na_balance", ylim = c(-0.1, 0.1), lwd=2, ylab = "meq/min", xlab = "Day")
plot(x4[,"time"]/24,x4[,"Na_concentration"], type="l") #, main = "Na_balance", ylim = c(-0.1, 0.1), lwd=2, ylab = "meq/min", xlab = "Day")

plot(x4[,"time"]/24,x4[,"plasma_renin_activity"], type="l")
plot(x4[,"time"]/24,x4[,"postglomerular_pressure"], type="l")
plot(x4[,"time"]/24, x4[,"GFR_ml_min"], type="l")

plot(x4[,"time"]/24, x4[,"arterial_pressure"], type="l")
plot(x4[,"time"]/24, x4[,"blood_volume_L"], type="l")
plot(x4[,"time"]/24, x4[,"e_cd_sodreab"], type="l")
