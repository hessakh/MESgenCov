
#' @keywords internal
appendPre <- function(site,preCSVs,obs,obsi,siObs,bi){
	if(bi == 1){
		site[,5+obsi] <- (10^-(site[,4+obsi]))*100794 #convert pH to H+ mg/L
		colnames(site)[5 + siObs] <- "H+ Conc"
	} 
		for (i in 1:(dim(site)[1])){
      		preTemp <- preCSVs[preCSVs$starttime>=site[i,2],]
      		preTemp <- preTemp[preTemp$endtime<=site[i,3],]
      		p       <- sum(preTemp$amount) 
      	site[i,5 + siObs+bi] <- (p)             #record weekly preciptation amt
      	#record chemical volume #change index here bases on obsi
      	site[i,6 + siObs+bi] <- ((p)*site[i,4+obsi+bi]) #get H+ Vol
    	}
    	colnames(site)[5 + siObs+bi] <- "preciptation"
    	if(bi){
    	  colnames(site)[6 + siObs+bi] <-"H+ Vol"
    	}else{colnames(site)[6 + siObs+bi] <- paste(obs[obsi],"Vol",sep="")}
    return(site)
}
