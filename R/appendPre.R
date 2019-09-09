
#' @keywords internal

appendPre <- function(site,preCSVs,obs,obsi,siObs,binary){
	j <- 1
    k <- 1
	if(binary == 1){
		for (i in 1:(dim(site)[1])){
    	#get out of smaller loop if current data is outside of current week range
      		if(preCSVs[j,3] > site[nrow(site),4]){
        		break}
      		p <- 0
      		for (j in k:(dim(preCSVs)[1])){
        		if (preCSVs[j,4] > site[i,4]){
          			k <- j
          			break
        		} else if(preCSVs[j,3] >= site[i,3]){
          		p <- ((p) + ((preCSVs[j,5])))
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
      		if(preCSVs[j,3] > site[nrow(site),4]){
        		break}
      		p <- 0
      		for (j in k:(dim(preCSVs)[1])){
        			if (preCSVs[j,4] > site[i,4]){
          		k <- j
          		break
        			} else if(preCSVs[j,3] >= site[i,3]){
          		p <- ((p) + ((preCSVs[j,5])))
        		}
      		}
      	site[i,6 + siObs] <- (p)             #record weekly preciptation amt
      	#record SO4 volume #change index here bases on obsi
      	site[i,7 + siObs] <- ((p)*site[i,5+obsi])
    	}
    	colnames(site)[6 + siObs] <- "preciptation"
    	colnames(site)[7 + siObs] <- paste(obs[obsi],"Vol",sep="")
	}
    site
}
