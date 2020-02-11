#'@keywords internal
#'This function takes a dataframe with weekly precipitation and concentration values,
#'
aggregateMonthly <-
function (binary,site,siObs,obs,obsi,totT,strtYrMo){#binary con or pH
   	  sitem <- data.frame(matrix(ncol=3, nrow=totT))
	    yrm        <- site[1,5]
	    strtYrMo   <- strtoi(strtYrMo)
	    sitem[1,1] <- strtYrMo
	    mStrt      <- ((strtYrMo) - floor(strtYrMo/100)*100)
	    yrdiff     <- (floor(yrm/100)- floor(strtYrMo/100))
	    mdiff      <- (yrm - floor(yrm/100)*100) - ((strtYrMo) - floor(strtYrMo/100)*100)
	    diffYrm    <- 12*yrdiff + mdiff

	    sitem[,2]  <- 1:totT   #t
      #fill in yrmonth
	    for(i in 1:totT){
	      sitem[i,1] <- strtYrMo + 100*(floor((i+mStrt-2)/12)) + ((i+mStrt-2)%%12)
	      sitem[i,3] <- NA
	    }
      #fill in deposit values
	    sitem[diffYrm + 1 ,3] <- 0   #SO4 vol
	    lyrm <- site[dim(site)[1],5] #last month

	    j = diffYrm + 1
	    c <- 1 #counter of weeks in month
	    ma <- site[1,7 + siObs] #monthly aggregate

    	if (binary == 1){
			for (i in 2:dim(site)[1]){
		      if(site[i,5] == yrm){
		        ma <- ma + site[i,7 + siObs]
		        c <- c +1
		        if (yrm == lyrm && i == nrow(site)){
		          sitem[j,3] <- ma/c
		        }
		      } else if(is.na(match(yrm,site[,5]))){
		        print(paste0("1yrm =  ",yrm))
		        yrm <- site[i,5]
		        print(paste0("1site =  ",site[i,5]))
		      }else if(site[i,5] > yrm){
		      	#monthly aggregate / recorded perciptation over the month
		        sitem[j,3] <- ma/c
		        c <- 1
		        ma <- site[i,7 + siObs]
		        yrm <- site[i,5]
		        j <- match(yrm, sitem[,1])
		      }
		    }
    	}else{
	    for (i in j:nrow(site)){
	      yrdiff     <- floor(site[i,5]/100)-floor(yrm/100)
	      if(site[i,5] == yrm){
	        ma <- ma + site[i,7+siObs]
	        c <- c +1
	        if (yrm == lyrm && i == nrow(site)){
	          sitem[j,3] <- ma/sum(site[(i-c+1):i,6 + siObs])
	        }
	      }else if(is.na(match(yrm,site[,5]))){
	        #print(paste0("1yrm =  ",yrm))
	        yrm <- site[i,5]
	        #print(paste0("1site =  ",site[i,5]))
	      } else if(site[i,5] > yrm){ #if row yrm is > curr yrm wrap up last yrm
	      	#monthly aggregate / recorded perciptation over the month
	      	if(sum(site[(i-c):(i-1),6 + siObs]) == 0){sump = NA}
	      	else{sump = sum(site[(i-c):(i-1),6 + siObs])} #sum of precepitation in month
	        sitem[j,3] <- ma/sump
	        #print(paste0("j =  ",j))
	        #print(paste0("sitem[j, 3] =  ",sitem[j,3]))
	        c <- 1
	        ma  <- site[i,7 + siObs]
	        yrm <- site[i,5]
	        j <- match(yrm, sitem[,1])
	      }
	    }
	}
	colnames(sitem) <- c("yrmonth", "t", paste(obs[obsi],"Con",sep=""))
	sitem
}
