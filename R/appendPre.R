
#' @keywords internal
appendPre <- function(site,preCSVs,obs,obsi,siObs,bi){
	j <- 1
  k <- 1
	if(bi == 1){
		site[,5+obsi] <- (10^-(site[,4+obsi]))*100794 #convert pH to H+ mg/L
		colnames(site)[5 + siObs] <- "H+ Conc"
	}
		for (i in 1:(dim(site)[1])){
    	#get out of smaller loop if current data is outside of current week range
      		if(preCSVs[j,2] > site[nrow(site),2]){ #col 2: startime, col 3: endtime
        		break}
      		p <- 0
      		for (j in k:(dim(preCSVs)[1])){
        		if (preCSVs[j,3] > site[i,3]){
          		k <- j
          		break
        		}else if(preCSVs[j,2] >= site[i,2]){
          		p <- p + preCSVs[j,4]
        		}
      		}#done looking for weekly p
      	site[i,5 + siObs+bi] <- (p)             #record weekly preciptation amt
      	#record SO4 volume #change index here bases on obsi
      	site[i,6 + siObs+bi] <- ((p)*site[i,4+obsi+bi]) #get H+ Vol
    	}
    	colnames(site)[5 + siObs+bi] <- "preciptation"
    	if(bi){
    	  colnames(site)[6 + siObs+bi] <-"H+ Vol"
    	}else{colnames(site)[6 + siObs+bi] <- paste(obs[obsi],"Vol",sep="")}
    return(site)
}
