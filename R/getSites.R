#' returns a vector of siteID's of sites with the most data for the specified time period for specific compunds or pH
#' @import lubridate
#' @import   stats
#' @importFrom utils data
#' @param startdateStr String, when to start analyzing data, format = "m/d/y H:M"
#' @param enddateStr   String, when to stop analyzing data, format = "m/d/y H:M"
#' @param maxn         Maximum number of sites required
#' @param mins         Minimum number of data required for each site
#' @param comp         Vector of strings, strings of compounds
#' @param optR         String, option of region "N","S","W",""
#' @return             Returns list of sites that have the largest amount of data for the given dates and compounds
#' @export
#' @examples getSites("01/01/83 00:00","12/31/86 00:00",30,102,"SO4","N")

getSites <- function(startdateStr,enddateStr,maxn,mins,comp,optR){
  #get data if it's not in the working directory
  if(!exists("weeklyConc") || !exists("preDaily")){
    try({utils::data("weeklyConc"); utils::data("preDaily")})
  }
  #check, still doesn't exist?
  if(!exists("weeklyConc") || !exists("preDaily")){
    message("Missing files),#running code that downloads necessary files from the NADP site")
    #getDataOffSite()
  }
  #get data if it's not in the working directory
  if(!exists("NADPgeo")){
    try({utils::data("NADPgeo",envir = environment());load("NADPgeo.rda",envir = environment())})
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
  rm(preCSV)
  
  obs    <- comp

  for (i in 1:length(obs)){
    obsiCSV <- match(obs[i],colnames(MESgenCov::weeklyConc))
    if(is.na(obsiCSV)){
      str1 <- paste(colnames(MESgenCov::weeklyConc[,5:14]), collapse = ", ")
      stop("Argument used in column comp of input data frame is not available. These are the options: ph, ", str1,".",collapse = " ")
    }
    cn <- colnames(conCSVf)
    conCSVf[,4+i] <- MESgenCov::weeklyConc[,obsiCSV]
    colnames(conCSVf) <- c(cn, obs[i])
    conCSVf <- conCSVf[conCSVf[,4+i]>-0.0001,] #filter out -/ive values
  }
  # #0.103 secs

  startdate <- as.POSIXct(startdateStr, format = "%m/%d/%y %H:%M")
  enddate   <- as.POSIXct(enddateStr  , format = "%m/%d/%y %H:%M")

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
  if(optR == "W"){
    siteDataCount <- siteDataCount[siteDataCount$long<=-100,]
  }else if(optR == "S"){
    siteDataCount <- siteDataCount[siteDataCount$long>=-100,]
    siteDataCount <- siteDataCount[siteDataCount$lat<=39.5,]
  }else if(optR == "N"){
    siteDataCount <- siteDataCount[siteDataCount$long>=-100,]
    siteDataCount <- siteDataCount[siteDataCount$lat>=39.5,]
  }
  
  #final list of sites
  if((dim(siteDataCount)[1]) < maxn && (dim(siteDataCount)[1])>0){
    finalSites <- siteDataCount$siteID
    my_list <- list("finalList" = finalSites, "data" = siteDataCount)
    message(paste0("Number of sites with data for inputted dates is less than specified max, number of sites found ",dim(siteDataCount)[1]))
    return(my_list)
    return(siteDataCount$siteID[1:(dim(siteDataCount)[1])])
  }else if ((dim(siteDataCount)[1]) < maxn) {
    message(paste0("Number of sites with data for inputted dates is less than specified max, number of sites found ",dim(siteDataCount)[1]))
    }else{
      finalSites <- siteDataCount$siteID[1:maxn]
      my_list <- list("finalList" = finalSites, "data" = siteDataCount, "startDate" = startdateStr,
                      "endDate" = enddateStr, "comp" = comp)
      return(my_list)
  }
}
