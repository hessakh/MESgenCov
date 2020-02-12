
#'gets weekly values and concentrations from daily precipitayion data
#'
#'@param binary Is pH being analyzed on the current loop? T/F (generated automatically in loop)
#'@param site   Data with compund/pH weekly values
#'@param siObs  Number of compunds analyzed
#'@param obs
#'@param obsi   Index of current compund being analyzed
#'@param startdate Startdate of analysis
#'@return sitew data that has weekly concentrations of compound, in the case of pH it gives exactly the weekly value plus a column with weekly precipitation values
#'@keywords internal

#@examples
# weeklyConc(F,site,1,1,as.POSIXct("01/01/83 00:00", format = "%m/%d/%y %H:%M"))

weeklyConc <-
  function (binary,site,siObs,obs,obsi,startdate){
   sitew <- data.frame(matrix(ncol=3, nrow=0))

	    sd <- site[1,3]
	    sitew[1,1] <- sd
	    class(sitew[1,1]) <- c('POSIXt','POSIXct')
	    currt <- round(as.numeric(difftime(as.POSIXct(sd,origin = origin), startdate, units = "weeks"))) + 1#current t,starts at 1
	    sitew[1,2] <- currt  #t # need to make sure this coincides with 1st week
	    sitew[1,3] <- 0          #comp vol
	    lsd <- round(as.numeric(difftime(as.POSIXct(site[dim(site)[1],3],origin = origin), startdate, units = "weeks"))) + 1#last week start date
	    j  <- 1 #tracks rows in sitew
	    c  <- 1 #counter of weeks in month
	    wa <- site[1,7 + siObs] # aggregate of compund

    	if (binary == 1){
			for (i in 2:dim(site)[1]){
				nextt <- round(as.numeric(difftime(as.POSIXct(site[i,3],origin = origin), startdate, units = "weeks")))+1#next t
		      if(nextt == currt){
		        wa <- wa + site[i,7 + siObs]
		        c  <- c + 1
		        if (currt == lsd && i == nrow(site)){
		          sitew[j,3] <- wa/c
		        }
		      } else if(nextt > currt){
		      	#weekly aggregate / recorded perciptation over the month
		        sitew[j,3] <- wa/c
		        c  <- 1
		        wa <- site[i,7 + siObs]
		        j  <- j+1
		        sitew[j,1] <- site[i,3]
		        #store t
		        currt <- nextt
				sitew[j,2] <- currt
		      }
		    }
    	}else{
	    for (i in 2:dim(site)[1]){
	    	nextt <- round(as.numeric(difftime(as.POSIXct(site[i,3],origin = origin), startdate, units = "weeks")))+1#next t
	      if(nextt == currt){
	        wa <- wa + site[i,7+siObs]
	        c  <- c + 1
	        if (currt == lsd && i == nrow(site)){
	          sitew[j,3] <- wa/sum(site[(i-c+1):i,6 + siObs])
	        }
	      } else if(nextt > currt){
	      	#weekly aggregate / recorded perciptation over the week
	      	if(sum(site[(i-c):(i-1),6 + siObs]) == 0){sump = 1}
	      	else{sump = sum(site[(i-c):(i-1),6 + siObs])}
	        sitew[j,3] <- wa/sump
	        c <- 1
	        wa <- site[i,7+siObs]
	        currt <- nextt
	        j = j+1
	        sitew[j,1] <- as.POSIXct(site[i,3],origin = origin)
	        #store t
	        sitew[j,2] <- currt
	      }
	    }
	}
	colnames(sitew) <- c("week starting" , "t", paste(obs[obsi],"Con",sep=""))
	sitew
}


