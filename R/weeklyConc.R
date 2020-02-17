#'gets weekly values and concentrations from daily precipitayion data
#'@keywords internal
#'
#'@param bi Is pH being analyzed on the current loop? T/F (generated automatically in loop)
#'@param site   Data with compund/pH weekly values
#'@param siObs  Number of compunds analyzed
#'@param obs    Vectpr of comp
#'@param obsi   Index of current compund being analyzed
#'@param startdate Startdate of analysis
#'
#'@return sitew data that has weekly concentrations of compound, in the case of pH it gives exactly the weekly value plus a column with weekly precipitation values

weeklyConc <-
  function (bi,site,siObs,obs,obsi,startdate){
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
	    wa <- site[1,7 + siObs+bi] # aggregate of compund

	    for (i in 2:dim(site)[1]){
	      nextt <- round(as.numeric(difftime(as.POSIXct(site[i,3],origin = origin), startdate, units = "weeks")))+1#next t
	      if(nextt == currt){
	        wa <- wa + site[i,7+siObs+bi]
	        c  <- c + 1
	        if (currt == lsd && i == nrow(site)){
	          sitew[j,3] <- wa/sum(site[(i-c+1):i,6 + siObs+bi])
	        }
	      } else if(nextt > currt){
	        #weekly aggregate / recorded perciptation over the week
	        if(sum(site[(i-c):(i-1),6 + siObs+bi]) == 0){sump = 1}
	        else{sump = sum(site[(i-c):(i-1),6 + siObs+bi])}
	        sitew[j,3] <- wa/sump
	        c <- 1
	        wa <- site[i,7+siObs+bi]
	        currt <- nextt
	        j = j+1
	        sitew[j,1] <- as.POSIXct(site[i,3],origin = origin)
	        #store t
	        sitew[j,2] <- currt
	      }
	    }
	colnames(sitew) <- c("week starting" , "t", paste(obs[obsi],"Con",sep=""))
	sitew
}


