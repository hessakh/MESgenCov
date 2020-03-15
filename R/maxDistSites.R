## get data

#' returns a vector of siteID's of sites with the most data for the specified time period for specific compunds or pH
#' @import lubridate
#' @import   stats
#' @importFrom utils data
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
  #get data if it's not in the working directory
  if(!exists("weeklyConc") || !exists("preDaily")){
    try({utils::data("weeklyConc",envir = environment()); utils::data("preDaily",envir = environment())})
  }
  #check, still doesn't exist?
  if(!exists("weeklyConc") || !exists("preDaily")){
    message("Missing files")#, running code that downloads necessary files from the NADP site")
    #getDataOffSite()
  }
  #get data if it's not in the working directory
  if(!exists("NADPgeo")){
    try({utils::data("NADPgeo",envir = environment())})
  }
  #check, still doesn't exist?
  if(!exists("NADPgeo")){
    message("Location data missing")
  }
  
  preCSV <- stats::na.omit(MESgenCov::preDaily)
  geoCSV <- MESgenCov::NADPgeo[,-2]
  colnames(geoCSV) <- c("siteID","lat","long")
  
  #filter out unnecessary columns
  conCSVf  <- MESgenCov::weeklyConc[,-5:-14]
  conCSVf  <- stats::na.omit(conCSVf)
  preCSVf  <- preCSV[preCSV$amount>-0.0001,]   #filter out -/ive values
  rm(preCSV) # for efficiency
  
  obs    <- comp
  
  #add back desired columns based on input in obs e.g add back SO4, pH, NO3 etc.
  ###Note that this will only keep data for dates where all comp data is present
  ###i.e. if comp = c("SO4","NO3") then only dates for which both elements have data will be included
  #tic()
  for (i in 1:length(obs)){
    obsiCSV <- match(obs[i],colnames(MESgenCov::weeklyConc))
    if(is.na(obsiCSV)){
      str1 <- paste(colnames(MESgenCov::weeklyConc[,5:14]), collapse = ", ")
      stop("Argument used in column comp of input data frame is not available. These are the options: ph, ", str1,".", collapse = " ")
    }
    cn <- colnames(conCSVf)
    conCSVf[,4+i] <- MESgenCov::weeklyConc[,obsiCSV]
    colnames(conCSVf) <- c(cn, obs[i])
    conCSVf <- conCSVf[conCSVf[,4+i]>-0.0001,] #filter out -/ive values
  }
  
  #format input
  startdate <- as.POSIXct(startdateStr, format = "%m/%d/%y %H:%M")
  enddate   <- as.POSIXct(enddateStr  , format = "%m/%d/%y %H:%M")
  
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
  
  #get sites from each filtered data set
  conCSVf$siteID <- toupper(conCSVf$siteID)
  preCSVf$siteID <- toupper(preCSVf$siteID)
  sitesCon <- unique(conCSVf$siteID)
  sitesPre <- unique(preCSVf$siteID)
  
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
  finalList     <- siteDataCount[startingSite,] #adds site with the most data
  remainSDC     <- siteDataCount[-startingSite,]
  trueMax       <- min(dim(siteDataCount)[1],maxn)
  
  if(trueMax < maxn){
    stop(paste("Number of sites with data for inputted dates is less than specified max, number of sites found",
               dim(siteDataCount)[1]))}
#make distance matrix
  distMat <- data.frame(matrix(ncol=length(siteDataCount$siteID), nrow=length(siteDataCount$siteID)))
  rownames(distMat) <- siteDataCount$siteID
  colnames(distMat) <- siteDataCount$siteID
  for(g in 1:length(siteDataCount$siteID)){
    for(h in g:length(siteDataCount$siteID)){
      if(g == h){
        distMat[g,g] <- 10^7
      }else{
        distMat[g,h] <- sqrt((siteDataCount[g,3]-siteDataCount[h,3])^2 + (siteDataCount[g,4]-siteDataCount[h,4])^2)
        distMat[h,g] <-  distMat[g,h]
      }
    }
  }
#
  for(i in 1:maxn-1){
    tempMat <- distMat[finalList$siteID,remainSDC$siteID] #only consider dists to sites in list
    tempVec <- NULL
    for(j in 1:length(remainSDC$siteID)){
      #store min dist to finalSites in tempVec
      tempVec <- c(tempVec,min(tempMat[,j]))
    }
    keepSitei <- match(max(tempVec), tempVec)
    finalList <- rbind(finalList,remainSDC[keepSitei,])#ignore indicies here for now
    remainSDC <- remainSDC[-keepSitei,]
  }
  #code is done, add colum with maxmin dist and test per favore
  reList <- list("finalList" = finalList$siteID[1:maxn], "data" = siteDataCount, "startDate" = startdateStr,
                 "endDate" = enddateStr, "comp" = comp)
  return(reList)
}
