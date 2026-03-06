calcEndTime <- function(heart_rate, endtime){


beat_duration = 60/heart_rate #min/beat

finaltime = floor(endtime/beat_duration)*beat_duration

return(finaltime)
}