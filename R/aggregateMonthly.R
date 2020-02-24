#'This function takes a dataframe with weekly precipitation and concentration values,
#'
#'@keywords internal

aggregateMonthly <-
function (bi,site,siObs,obs,obsi,totT,strtYrMo){#binary con or pH
  sitem <- data.frame(matrix(ncol=3, nrow=totT))
  yrm        <- site[1,5]
  strtYrMo   <- strtoi(strtYrMo)
  mStrt      <- (strtYrMo)%%12
  
  #fill in t
  sitem[,2]  <- 1:totT
  #fill in yrmonth
  for(i in 1:totT){
    sitem[i,1] <- strtYrMo + 100*(floor((i+mStrt-2)/12)) + ((i+mStrt-2)%%12)
    sitem[i,3] <- NA
  }
  for (i in 0:(totT-2+mStrt)){# handles case with no precip data for month and skips in months
    siteyrmonth <- site[site$yrmonth==strtYrMo-mStrt+1+(i%%12)+((i)%/%12)*100,]
    if(sum(siteyrmonth[,obsi + 5+bi]) == 0){
      sitem[i+1,3] <- NA
    }else{
      sitem[i+1,3] <- (sum(siteyrmonth[,obsi+6+bi])/sum(siteyrmonth[,obsi+5+bi]))
    } #col obsi+7 is volume of comp, col obsi+6 is precipitation amount fone one comp analysis
  }
  if(bi){
    colnames(sitem) <- c("yrmonth", "t", paste("H+ ","Con",sep=""))
  }else{colnames(sitem) <- c("yrmonth", "t", paste(obs[obsi],"Con",sep=""))}
  sitem
}


