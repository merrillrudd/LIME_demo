## calculates the amount of time the gear was in the water

calc_soak_time <- function(time1, time2, hr=TRUE){
	if(time1==""|time2=="") return(NA)
	hr1 <- tryCatch(as.numeric(strsplit(time1, ":")[[1]][1]), error=function(e) NA)
	min1 <- tryCatch(as.numeric(strsplit(time1, ":")[[1]][2]), error=function(e) NA)

	hr2 <- tryCatch(as.numeric(strsplit(time2, ":")[[1]][1]), error=function(e) NA)
	min2 <- tryCatch(as.numeric(strsplit(time2, ":")[[1]][1]), error=function(e) NA)
	
	if(is.na(hr2)) return(NA)
	if(is.na(h1)) return(NA)
	if(hr2 < 12) hr2_v2 <- 24 + hr2
	if(hr2 >= 12) hr2_v2 <- hr2

	hr_elapsed <- hr2_v2 - hr1
	min_elapsed <- min2 - min1

	min_return <- hr_elapsed*60 + min_elapsed
	hr_return <- round(min_return/60,2)

	if(hr==TRUE) return(hr_return)
	if(hr==FALSE) return(min_return)
}