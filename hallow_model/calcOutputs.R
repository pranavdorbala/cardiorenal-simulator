calcOutputs <- function(x){
  
  library(plyr)
  #source("calcByPeriod.R")
  
  #Calculate LV ESV, SV, and EV
  x=data.frame(x)
  
  x = calcByPeriod(x, "min", "LV_volume", "LV_ESV" )

  x$stroke_volume_ml = (x$LV_EDV - x$LV_ESV)*1e6
  x$ejection_fraction = (x$LV_EDV - x$LV_ESV)/x$LV_EDV
  
  
  x = calcByPeriod(x, "min", "arterial_pressure", "DBP" )
  x = calcByPeriod(x, "max", "arterial_pressure", "SBP" )
  
  #plot(x[,"arterial_pressure"], type="l")
  #points(x[,"SBP"], type="l", col="red")
  #points(x[,"DBP"], type="l", col="red")
  
  x$SBP_mmHg = x$SBP*0.0075
  x$DBP_mmHg = x$DBP*0.0075
  
  x = calcByPeriod(x,"max", "LV_pressure", "LV_PSP")
  x = calcByPeriod(x,"total", "mitral_valve_leak", "MV_leak")
  
 return(x) 
}


calcByPeriod <- function(x, functype, varname, newvarname) {
  

 
  
  if (functype == "total") {
   
    thisdata = x[, names(x) %in% c("periods",varname, "beat_time")]
    thisvar = ddply(thisdata, "periods", function(thisPeriod) {
      lastVal = thisPeriod[,varname][thisPeriod$beat_time == min(thisPeriod$beat_time)]
      periodVal = thisPeriod[,varname] - lastVal
      data.frame(dv = periodVal)
    })
    names(thisvar)[names(thisvar) == "dv"] =  newvarname
    x[,newvarname] =   thisvar[,newvarname]
  } else {
    
    thisdata = x[, names(x) %in% c("periods",varname)]
    names(thisdata)[names(thisdata) != "periods"] = "value"
    
    if (functype == "max") {
      thisvar = ddply(thisdata, "periods", summarise, dv = max(value))
    } else if (functype == "min") {
      thisvar = ddply(thisdata, "periods", summarise, dv = min(value))
    }
    names(thisvar)[names(thisvar) == "dv"] =  newvarname
    x1 = merge(x,thisvar, by="periods")
    x[,newvarname] =   x1[,newvarname]
    }

  
  


  #take care of the first beat - if a partial beat, use value from next beat
  if (x$beat_time[x$periods == min(x$periods)] > 0.1) {
    x[x$periods == min(x$periods),newvarname] = x[x$periods == min(x$periods) + 1, newvarname][1]
  }
  
    #take care of the last beat - if a partial beat, use value from previous beat
  if (x$beat_time[x$periods == max(x$periods)] < 0.9) {
    x[x$periods == max(x$periods),newvarname] = x[x$periods == max(x$periods) - 1, newvarname][1]
  }
  
  return(x)
}


