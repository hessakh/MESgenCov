
#' @keywords internal

appendPre <- function(site,preCSVs,obs,obsi,siObs,binary){
	j <- 1
    k <- 1
	if(binary == 1){
		for (i in 1:(dim(site)[1])){
    	#get out of smaller loop if current data is outside of current week range
      		if(preCSVs[j,3] > site[nrow(site),4]){#j here should be k but test to find out
        		break}
      		p <- 0
      		for (j in k:(dim(preCSVs)[1])){
        		if (preCSVs[j,4] > site[i,4]){ #endtime (day) is after endtime for weekly data
          			k <- j                     #store index in preCSV to start from for the next row of data
          			break
        		} else if(preCSVs[j,3] >= site[i,3]){ #starttime is within starttime of weekly data
          		p <- ((p) + ((preCSVs[j,5])))       # add precipitation value to current store of weekly precipitation
        		}
      		}
      		site[i,6 + siObs] <- (p)             #record weekly preciptation amt
      		#record SO4 volume #change index here bases on obsi
      		site[i,7 + siObs] <- site[i,5+obsi]
      	}
      	colnames(site)[6 + siObs] <- "preciptation"
    	colnames(site)[7 + siObs] <- obs[obsi]
	}else{
		for (i in 1:(dim(site)[1])){
    	#get out of smaller loop if current data is outside of current week range
      		if(preCSVs[j,3] > site[nrow(site),4]){ #col 3: startime, col 4: endtime
        		break}
      		p <- 0
      		for (j in k:(dim(preCSVs)[1])){
        		if (preCSVs[j,4] > site[i,4]){
          		k <- j
          		break
        		}else if(preCSVs[j,3] >= site[i,3]){
          		p <- ((p) + ((preCSVs[j,5])))
        		}
      		}#done looking for weekly p
      	site[i,6 + siObs] <- (p)             #record weekly preciptation amt
      	#record SO4 volume #change index here bases on obsi
      	site[i,7 + siObs] <- ((p)*site[i,5+obsi])
    	}
    	colnames(site)[6 + siObs] <- "preciptation"
    	colnames(site)[7 + siObs] <- paste(obs[obsi],"Vol",sep="")
	}
    return(site)
}
