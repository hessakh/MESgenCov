## get data

#' returns a vector of siteID's of sites with the most data for the specified time period for specific compunds or pH
#' @import lubridate
#' @param startdateStr String, when to start analyzing data, format = "m/d/y H:M"
#' @param enddateStr   String, when to stop analyzing data, format = "m/d/y H:M"
#' @param maxn         Number of sites required
#' @param mins         Minimum number of data required for each site
#' @param comp         Vector of strings, strings of compounds
#' @param startingSite Integer, largest sample index of initial site to include
#' @return             Returns list of sites that have the largest amount of data for the given dates and pollutants
#' @export
#' @examples maxDistSites("01/01/83 00:00","12/31/86 00:00",50,100,"SO4",1)

maxDistSites <- function(startdateStr,enddateStr,maxn,mins,comp,startingSite){
  #tic()
  #get data if it's not in the working directory
  if(!exists("weeklyCSV") || !exists("preDailyCSV")){
    try({data("weeklyCSV"); data("preDailyCSV")})
  }
  #check, still doesn't exist?
  if(!exists("weeklyCSV") || !exists("preDailyCSV")){
    message("Missing files, running code that downloads necessary files from the NADP site")
    getDataOffSite()
  }

  #get data if it's not in the working directory
  if(!exists("NADPgeo")){
    try({data("NADPgeo")})
  }
  #check, still doesn't exist?
  if(!exists("NADPgeo")){
    message("Location data missing")
  }

  conCSV <- weeklyCSV
  preCSV <- na.omit(preDailyCSV)
  geoCSV <- NADPgeo
  colnames(geoCSV) <- c("siteID","city","lat","long")
  #toc() #0.914 sec if data is already loaded if not 26 sec
  #tic()
  #filter out unnecessary columns
  conCSVf  <- conCSV[,-6:-31]
  conCSVf  <- na.omit(conCSVf)
  preCSVf  <- preCSV[,-1]
  preCSVf  <- preCSVf[preCSVf$amount>-0.0001,]   #filter out -/ive values

  obs    <- comp
  #toc() #0.65 secs
  #add back desired columns based on input in obs e.g add back SO4, pH, NO3 etc.
  ###Note that this will only keep data for dates where all comp data is present
  ###i.e. if comp = c("SO4","NO3") then only dates for which both elements have data will be included
  #tic()
  for (i in 1:length(obs)){
    obsiCSV <- match(obs[i],colnames(conCSV))
    cn <- colnames(conCSVf)
    conCSVf[,5+i] <- conCSV[,obsiCSV]
    colnames(conCSVf) <- c(cn, obs[i])
    conCSVf <- conCSVf[conCSVf[,5+i]>-0.0001,] #filter out -/ive values
  }
  #toc() #0.103 secs
  #tic()
  startdate <- as.POSIXct(startdateStr, format = "%m/%d/%y %H:%M")
  enddate   <- as.POSIXct(enddateStr  , format = "%m/%d/%y %H:%M")
  #toc()
  #transform date into date format w/ time

  #filter by date
  strtYrMo  <- format(startdate,"%Y%m")
  endYrMo   <- format(enddate,  "%Y%m")
  conCSVf   <- conCSVf[conCSVf$yrmonth>=strtYrMo,]
  conCSVf   <- conCSVf[conCSVf$yrmonth<=endYrMo,]

  #filter dates for precipitation data
  d1 <- startdate
  d2 <- enddate
  preCSVf <- preCSVf[preCSVf$starttime >=d1,]
  preCSVf <- preCSVf[preCSVf$endtime   <=d2,]

  #get sites from each data set
  conCSV$siteID <- toupper(conCSV$siteID)
  preCSV$siteID <- toupper(preCSV$siteID)
  sitesCon <- unique(conCSV$siteID)
  sitesPre <- unique(preCSV$siteID)

  #get common sites from weekly and daily data
  commonSites <- intersect(sitesCon,sitesPre)
  # not included sites
  setDiff <- union(setdiff(sitesCon,sitesPre),setdiff(sitesPre,sitesCon))

  #count data for each site
  siteDataCount <- data.frame(matrix(ncol=2, nrow=0))
  k <- 0 #skip counter
  for(i in seq_along(commonSites)){
    tempCSV <- conCSVf[conCSVf$siteID == toString(commonSites[i]),]
    if ((dim(tempCSV)[1]) >= mins){
      siteDataCount[i-k,1] <- commonSites[i]
      siteDataCount[i-k,2] <- (dim(tempCSV)[1])
    } else{
      k <- k + 1
    }
  }
  colnames(siteDataCount) <- c("siteID","count")
  siteDataCount <- merge(siteDataCount, geoCSV, by="siteID")
  siteDataCount <- siteDataCount[order(siteDataCount$count,decreasing = TRUE),]
  finalList <- siteDataCount[startingSite,]
  remainSDC <- siteDataCount[-startingSite,]
  trueMax <- min(dim(siteDataCount)[1],maxn)
  #print(dim(siteDataCount)[1])
  if(trueMax < maxn){#toc()
    stop(paste("Number of sites with data for inputted dates is less than specified max, number of sites found",
                  dim(siteDataCount)[1]))}
  for(i in 1:trueMax-1){
    l1 <- dim(remainSDC)[1]
    for(j in 1:l1){ #one more for loop here to sum through final list, maybe switch up loops
      l2 <- dim(finalList)[1]
      partialDist <- 0
      tempVec <- NULL
      for(m in 1:l2){
        tempVec <- c(tempVec,sqrt((finalList[m,4]-remainSDC[j,4])^2 + (finalList[m,5]-remainSDC[j,5])^2))
        remainSDC[j,6] <- partialDist  + sqrt((finalList[m,4]-remainSDC[j,4])^2 + (finalList[m,5]-remainSDC[j,5])^2)
        partialDist    <- remainSDC[j,6]
      }
      remainSDC[j,6] <- min(tempVec)
      #print(remainSDC[1:5,])
    }
    colnames(remainSDC) <- c("siteID", "count", "city", "lat", "long", "dist")
    remainSDC <- remainSDC[order(remainSDC$dist,decreasing = TRUE),]
    finalList <- rbind(finalList,remainSDC[1,1:5])
    remainSDC <- remainSDC[-1,]
  }
  finalList <- finalList[1:maxn,]
  reList <- list("finalList" = finalList$siteID, "data" = siteDataCount, "startDate" = startdateStr,
                 "endDate" = enddateStr, "comp" = comp)
  return(reList)
}

