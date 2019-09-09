aggregateMonthly <-
function (binary,site,siObs,obs,obsi){#binary con or pH
   	    sitem <- data.frame(matrix(ncol=0, nrow=0))

	    yrm <- site[1,5]
	    sitem[1,1] <- yrm
	    sitem[1,2] <- 1   #t
	    sitem[1,3] <- 0   #SO4 vol
	    lyrm <- site[dim(site)[1],5] #last month

	    j = 1
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
		      } else if(site[i,5] > yrm){
		      	#monthly aggregate / recorded perciptation over the month
		        sitem[j,3] <- ma/c
		        c <- 1
		        ma <- site[i,7 + siObs]
		        yrm <- site[i,5]
		        j = j+1
		        sitem[j,1] <- yrm
		        #store t
		        if(sitem[j,1] - sitem[j-1,1] >= 89){#deals with jumps in year
		          sitem[j,2] <- sitem[j,1] - sitem[j-1,1] + sitem[j-1,2] -89 +1
		        }else{
		          sitem[j,2] <- sitem[j,1] - sitem[j-1,1] + sitem[j-1,2] #tracks month
		        }
		      }
		    }
    	}else{
	    for (i in 2:dim(site)[1]){
	      if(site[i,5] == yrm){
	        ma <- ma + site[i,7+siObs]
	        c <- c +1
	        if (yrm == lyrm && i == nrow(site)){
	          sitem[j,3] <- ma/sum(site[(i-c+1):i,6 + siObs])
	        }
	      } else if(site[i,5] > yrm){
	      	#monthly aggregate / recorded perciptation over the month
	      	if(sum(site[(i-c):(i-1),6 + siObs]) == 0){sump = 1}
	      	else{sump = sum(site[(i-c):(i-1),6 + siObs])}
	        sitem[j,3] <- ma/sump
	        c <- 1
	        ma <- site[i,7+siObs]
	        yrm <- site[i,5]
	        j = j+1
	        sitem[j,1] <- yrm
	        #store t
	        if(sitem[j,1] - sitem[j-1,1] >= 89){#deals with jumps in year
	          sitem[j,2] <- sitem[j,1] - sitem[j-1,1] + sitem[j-1,2] -89 +1
	        }else{
	          sitem[j,2] <- sitem[j,1] - sitem[j-1,1] + sitem[j-1,2] #tracks month
	        }
	      }
	    }
	}
	colnames(sitem) <- c("yrmonth", "t", paste(obs[obsi],"Con",sep=""))
	sitem
}
