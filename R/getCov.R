#'Creates covariance matrix from normalized NADP monitor data
#'gamma version with optional inputs
#'@import   MVN
#'@importFrom rmatio write.mat
#'@importFrom EnvStats rosnerTest
#'@importFrom utils data
#'@param df data frame of input arguments specified in vignette and default set of inputs are in data("defultInput")

#'@return   List of model summaries at each site, covariance matrix, sites analyzed, predicted values by univariate models, data frame of residuals and more (see vignette for full list)
#'@export
#'@examples
#'   getCov(structure(list(
#' startdateStr = "01/01/83 00:00", enddateStr   = "12/31/86 00:00", 
#'comp = "SO4", use36 = TRUE, siteAdd = NULL, outlierDatesbySite = NULL,
#'siteOutliers = NULL,  removeOutliers = NULL, plotMulti = FALSE,  sitePlot = NULL,
#'plotAll = FALSE, writeMat = FALSE, seas = 12, r = 1, k = 1),
#'.Names = c("startdateStr","enddateStr","comp","use36","siteAdd",
#'"outlierDatesbySite","siteOutliers","removeOutliers","plotMulti","sitePlot",
#'"plotAll","writeMat","seas","r","k")  ,row.names = c(NA, -1L) 
#',class = "data.frame"))


getCov <- function(df){
  #get data if it's not in the working directory
  if(!exists("weeklyConc") || !exists("preDaily")){
    try({utils::data("weeklyConc",envir = environment()); utils::data("preDaily",envir = environment())})
  }
  #check, still doesn't exist?
  if(!exists("weeklyConc") || !exists("preDaily")){
    stop("Missing files, running code that downloads necessary files from the NADP site")
  }
  #check for multiple pollutants
  if(length(df$comp) > 1){
    stop("Package can only model data from one pollutant/observed variable")
    #later versions of package will be able to handle this
  }
  
  ##code for optional input functionality
  dfInp <- df
  weeklyB      <- FALSE#df$weeklyB
  startdateStr <- df$startdateStr
  enddateStr   <- df$enddateStr
  comp         <- df$comp
  use36        <- df$use36
  siteAdd      <- df$siteAdd
  outlierDatesbySite <- df$outlierDatesbySite
  siteOutliers <- df$siteOutliers
  removeOutlier <-df$removeOutliers 
  plotMulti    <- df$plotMulti
  sitePlot     <- df$sitePlot
  plotAll      <- df$plotAll
  writeMat     <- df$writeMat
  seas         <- df$seas
  r            <- df$r
  k            <- df$k
  
  if (is.null(siteOutliers)){
    showOutliers = FALSE
  }else{showOutliers = TRUE}
  if (is.null(sitePlot)){
    plotB = FALSE
  }else{plotB = TRUE}
  
  ##### store data
  
  conCSV <- MESgenCov::weeklyConc
  preCSVf <- stats::na.omit(MESgenCov::preDaily)
  conCSV$siteID  <- toupper(conCSV$siteID) #combine data
  preCSVf$siteID <- toupper(preCSVf$siteID)
  
  #filter out unnecessary columns
  conCSVf  <- conCSV[,-5:-14]
  conCSVf  <- stats::na.omit(conCSVf)
  
  startdate  <- as.POSIXct(startdateStr, format = "%m/%d/%y %H:%M")
  enddate    <- as.POSIXct(enddateStr  , format = "%m/%d/%y %H:%M")
  strtYrMo   <- format(startdate,"%Y%m")
  endYrMo    <- format(enddate,  "%Y%m")
  diffYrm    <- mondf(startdate,enddate)+1
  nonneg <- 0
  if(weeklyB){
    totT <- as.integer(difftime(enddate,startdate,units="weeks"))
    nonneg <- 1
    if(is.null(seas)){seas = 52}
  }else if(!(floor(mondf(startdate,enddate))==(mondf(startdate,enddate)))){
    stop("Number of months is not integer")
  }else if(diffYrm <= 1){
    stop("Number of months is less than or equal to one")
  }else{
    totT <- (mondf(startdate,enddate))+1
    if(is.null(seas)){seas = 12}
  }
  
  #add default sites
  if(use36){
    cat36 <- c("AL10","IL11","IL18","IL19","IL35","IL47","IL63","IN34",
               "IN41","MA01","MA13","MD13","MI09","MI26","MI53",
               "NC03","NC34","NC41","NJ99","NY08","NY10","NY20",
               "NY52","NY65","OH17","OH49","OH71","PA15","PA29",
               "PA42","TN00","TN11","VA13","VT01","WI28","WV18")#from Guttorp, Le 1992
  }else{cat36 <- NULL}
  cati   <- union(cat36,siteAdd[[1]])
  obs    <- comp
  
  #add back desired columns based on input in obs e.g add back SO4, pH, NO3 etc.
  for (i in 1:length(obs)){
    obsiCSV           <- match(obs[i],colnames(conCSV))
    cn                <- colnames(conCSVf)
    conCSVf[,4+i]     <- conCSV[,obsiCSV] #
    colnames(conCSVf) <- c(cn, obs[i])
  }
  rm(conCSV)
  
  #filter by date
  conCSVf   <- conCSVf[conCSVf$yrmonth>=strtYrMo,]
  conCSVf   <- conCSVf[conCSVf$yrmonth<=endYrMo,]
  
  #filter dates for precipitation data
  d1 <- startdate
  d2 <- enddate
  preCSVf <- preCSVf[preCSVf$starttime >=d1,]
  preCSVf <- preCSVf[preCSVf$endtime   <=d2,]
  
  preCSVf$amount    <- as.numeric(preCSVf$amount) #change data type of column
  preCSVf <- preCSVf[preCSVf$amount>-0.0001,]     #filter out -/ive values

  
  #initialize
  
  rv    <- integer(5); if(r != 0){rv[1:r] <- rep(1L,r)}
  kv    <- integer(5); if(k != 0){kv[1:k] <- rep(1L,k)}
  kk    <- k
  
  tNA     <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp before na.omit
  tafNA   <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp
  y0      <- vector(mode="list", length=(length(cati)*length(obs))) #monthly/weekly concentraition values
  ylogNA  <- vector(mode="list", length=(length(cati)*length(obs))) #log data w/ NA values
  ylog    <- vector(mode="list", length=(length(cati)*length(obs))) #log data
  e    <- vector(mode="list", length=(length(cati)*length(obs))) #error
  e2   <- vector(mode="list", length=(length(cati)*length(obs))) #error with NA
  mods <- vector(mode="list", length=(length(cati)*length(obs))) #model summaries
  vpredl <- vector(mode="list", length=(length(cati)*length(obs))) #list of vector of predictions
  
  si   <- 1
  obsi <- 1
  siObs<- length(obs)
  oc   <- 1

  if(plotAll){
    graphics::par(mfrow = c(2,2),
        oma = c(0,0,0,0) + 0.1,
        mar = c(0,0,0,0) + 0.1)
  }
  for(s in 1:(length(cati)*length(obs))){
    #change obs
    if (s%%(length(cati)) == 1 && s != 1){obsi = obsi + 1; si = 1}
    
    #filter site (weekly concentration data)
    conCSVk <- conCSVf[conCSVf$siteID == toString(cati[si]),]
    preCSVk <- preCSVf[preCSVf$siteID == toString(cati[si]),]
    
    if(nrow(conCSVk)==0 || nrow(preCSVk)==0){
      #store parameters
      tNA[[s]]     <- NA
      to[[s]]     <- NA
      e[[s]]      <- NA
      e2[[s]]     <- NA
      ylog[[s]]   <- NA
      y0[[s]]     <- NA
      ylogNA[[s]] <- NA
      mods[[s]]   <- NA
      vpredl[[s]] <- NA
      message(paste0("Missing data for ", cati[si]," ",obs[obsi],
                     " check if site has data for inputted dates in data file weeklyConc, preDaily"))
    }else{
      #filter out missing data, -/ive values, and order dataframes by date
      site <- conCSVk
      site <- site[site[,4 + obsi]>-0.0001,]
      site <- site[order(site$dateon),]
      preCSVk <- preCSVk[order(preCSVk$starttime,decreasing=F),]

      
      #aaggregates weekly precipitation
      if (obs[obsi] == "ph"){bi = 1}else{bi = 0}
      site <- appendPre(site,preCSVk,obs,obsi,siObs,bi)
      site <- stats::na.omit(site)
      
      #aggregate precipitation data monthly and get concentration values
      if(!weeklyB){sitem <- aggregateMonthly(bi,site,siObs,obs,obsi,totT,strtYrMo)
      }#else{       sitem <- weeklyConcT(bi,site,siObs,obs,obsi,startdate)}
      
      if (length(outlierDatesbySite) != 0){
        if(outlierDatesbySite[oc] == cati[si]){
          for(i in (oc+1):length(outlierDatesbySite)){
            ifInt  <- suppressWarnings(as.integer(outlierDatesbySite[i]))
            if (!is.na(ifInt)){
              j     <- as.integer(outlierDatesbySite[i])
              index <- match(j, sitem$t)
              if (!is.na(index)){
                sitem[index,] <- NA
              }
            }else{break}
          }
          oc = i #outlier site counter #skips t
        }
      }
      p <- 1
      #deterministic univariate model for one site
      cn <- colnames(sitem)
      minsite3 <- min(stats::na.omit(sitem[,3]))
      if(!weeklyB || minsite3 >= 0.0001){
        minsite3 <- 0
        nonneg   <- 0
      }else if(minsite3 == 0){
        nonneg = 0.00001
      }
      sitem[,4] <- log(sitem[,3] + abs(minsite3) + nonneg)
      ylogNA <- sitem[,4]
      #moved y0, ylogNA from here to avoid NAs in plot call
      sitem     <- stats::na.omit(sitem)
      colnames(sitem) <- c(cn[1:2],"Conc" ,"log")
      y0[[s]]     <- sitem$Conc
      ylog[[s]]   <- sitem$log
      tafNA[[s]]  <- sitem$t #same as t for monthly different for weekly
      t <- sitem$t
      cyclicTrend <- (I(cos(t*(2*pi/seas))^p)   + I(sin(t*(2*pi/seas))^p))*kv[1]   +
                     (I(cos(t*(2*pi/seas)*2)^p) + I(sin(t*(2*pi/seas)*2)^p))*kv[2] +
                     (I(cos(t*(2*pi/seas)*3)^p) + I(sin(t*(2*pi/seas)*3)^p))*kv[3] +
                     (I(cos(t*(2*pi/seas)*4)^p) + I(sin(t*(2*pi/seas)*4)^p))*kv[4] +
                     (I(cos(t*(2*pi/seas)*5)^p) + I(sin(t*(2*pi/seas)*5)^p))*kv[5] +
                   t*rv[1] + (t^2)*rv[2] + (t^3)*rv[3] + (t^4)*rv[4] + (t^5)*rv[5]
      dfc  <- data.frame(cbind(sitem$log,cyclicTrend))
      dfc  <- data.frame(cbind(dfc,t))
      mod <- definelm(sitem$log,t,dfc,r,kk,seas)
      er  <- stats::residuals(mod)
      summary(mod)
      mods[[s]] <- summary(mod)
      tc       <- data.frame(1:totT)
      currResi <- data.frame(cbind(t,er))
      colnames(currResi) <- c("t","error")
      colnames(tc)      <- c("t")
      cR       <- merge(tc,currResi,by = "t", all = TRUE) #merge(dfRes,currRes, by = "t", all = F)
      e2[[s]]  <- cR[,2]
      
      #store predicted value vector
      options(warn = 2)
      if(length(tafNA[[s]]) < totT){
        vpred <- predict(mod,newdata = tc)
      }else{vpred <- predict(mod)}
      options(warn = 0)
      
      #fill in missing t ##there's a better way of doing this
      se1  <- sqrt(deviance(mod)/df.residual(mod))
      fi   <- 0 #to loop back
      maxfi<- totT - length(t)
      fe   <- 1:(maxfi) #fake end to keep loop going #ignore if t=48
      t    <- c(t,fe)
      y1   <- c(sitem$log,fe)
      er   <- c(er,fe)
      kl    <- 1

      for(i in 1:(totT+maxfi-1)){
        if(t[kl] != i-fi){ #if t is skipped
          cy1 <- vpred[i-fi]
          ry1 <- y1[kl:length(y1)]
          set.seed(i+s)
          e0 <- rnorm(1,0,se1)
          y1 <- c(y1[1:kl-1],cy1+e0,ry1)
          er <- c(er[1:kl-1],e0,er[kl:length(er)]) # change to stochastic #change to 0
          t  <- c(t[1:kl-1],i-fi,t[kl:length(t)])
          fi <- fi + 1
        }else{kl = kl + 1}
      }
      if(length(t)>totT){t<-t[1:totT]; y1 <- y1[1:totT];er<-er[1:totT]}
      
      #store parameters
      e[[s]]     <- unname(er)
      vpredl[[s]] <- vpred
      
      #plot all sites if indicated by user
      if(plotAll){
        if(s%%4 == 1 && s!=1){
          grDevices::dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
          graphics::par(mfrow = c(2,2))}
        graphics::par(mar=c(4,4,2,2))
        tc <- 1:length(vpred)
        graphics::plot(sitem$t,sitem$log, ylab="Log concentration",main = paste(cati[si],obs[obsi]), xlab = "t (months)")
        graphics::par(new=TRUE)
        graphics::lines(x = tc, y = vpred, col ="blue")
      }
    }
    si <- si + 1
  }
  ############################END BIG LOOP##############################
  
  options(warn=-1)
  dfRes  <- NULL
  dfRes2 <- NULL
  si    <- 2
  obsi  <- 1
  tc <- 1:totT
  #loop to create dataframe of all residuals
  dfRes <- cbind(tc,e[[1]])
  dfRes2 <- cbind(tc,e2[[1]]) #with NA values
  colnames(dfRes) <- c("t",paste(cati[1],obs[1], sep=""))
  colnames(dfRes2) <- c("t",paste(cati[1],obs[1], sep=""))
  for( i in 2:(length(cati)*length(obs))){
    if (i%%(length(cati)) == 1){obsi = obsi + 1; si = 1}
    if(!is.na(e[[i]])){
      currRes <- cbind(tc,e[[i]])
      colnames(currRes) <- c("t",paste0(cati[si],obs[obsi]))
      dfRes   <- merge(dfRes,currRes, by = "t", all = F)
    }
    currRes2 <- cbind(tc,e2[[i]])
    colnames(currRes2) <- c("t",paste0(cati[si],obs[obsi]))
    dfRes2   <- merge(dfRes2,currRes2, by = "t", all.x = TRUE,all.y = TRUE)
    si <- si+1
  }
  #create covariance matrix
  covxx <- data.frame(cov(dfRes[,-1]))
  options(warn=0)
  
  #produce multivariate analysis
  graphics::par(mar=c(1,1,1,1))
  MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "none", multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", showOutliers = FALSE,showNewData = FALSE)
  univariateTest <- MVDw$univariateNormality
  MVDw
  if(plotMulti){
    grDevices::dev.new(width = 4, height = 5, noRStudioGD = TRUE)
    #graphics::par(mfrow=c(1,2))
    MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25,
                alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW", 
                univariatePlot = "none", multivariatePlot = "qq", multivariateOutlierMethod = "none", 
                bc = FALSE, bcType = "rounded", showOutliers = FALSE, showNewData = FALSE)
    MVDw
  }
  #sitePlot
  if(plotB){
    tc <- 1:totT
    for (g in 1:length(sitePlot[[1]])){
      grDevices::dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
      graphics::par(mar=c(4,4,2,2))
      i = match(sitePlot[[1]][g], cati)
      if(!is.na(i)){
        graphics::plot(tafNA[[i]],ylog[[i]], ylab="Log concentration",main = toString(sitePlot[[1]][g]),xlab = "t (months)")
        graphics::par(new=TRUE)
        graphics::lines(x = tc, y = vpredl[[i]], col ="blue")
      }else{warning("Site in sitePlot was not in the vector of sites that were analyzed. Make sure the site ID in sitePlot is in the column siteAdd of the input data frame.")}
    }
  }
  #outlier analysis
  rosnerT   <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp
  siteRosner <- NULL
  noutliers  <- 0
  if(showOutliers == TRUE){
    if (length(siteOutliers[[1]])==1){
      #find outliers in each site
      #Store Site outliers
      i <- match(siteOutliers[[1]][1],cati)
      if(is.na(i)){
        message("siteOutliers string doesn't match any colnames, these are the options:")
        print(colnames(dfRes))
        siteRosner = NULL}
      else{siteRosner = rosnerTest(dfRes[,i+1])
          rosnerT[[i]] <- rosnerTest(dfRes[,i+1])}
    }else if (length(siteOutliers[[1]]) > 1){
      for (z in 1:length(siteOutliers[[1]])){
        i <- match(siteOutliers[[1]][z],cati)
        if(is.na(i)){
          message("siteOutliers string doesn't match any colnames, these are the options:")
          print(cati)
          siteRosner = NULL
        }else{siteRosner = c(siteRosner, rosnerTest(dfRes[,i+1]))
              rosnerT[[i]] <- rosnerTest(dfRes[,i+1])}
              noutliers    <- noutliers + rosnerTest(dfRes[,i+1])$n.outliers
      }
    }
  }
  if(!is.null(removeOutlier)){
    if(siteOutliers[[1]] %contain% removeOutlier[[1]]){
      if(length(removeOutlier[[1]])>=2){
        sitesOut <- removeOutlier[[1]]
        outlierDatesbySite <- takeOutAllOutliers(sitesOut,rosnerResult = rosnerT,cati)$outBySite
        outSites           <- takeOutAllOutliers(sitesOut,rosnerResult = rosnerT,cati)$sites
        if(length(outlierDatesbySite) >= 2){
          updatedPars        <- reEvaluateSites(dfInp, tafNA,startdate,enddate,totT,
                                                y0,ylog,e,e2,mods,vpredl, outlierDatesbySite,outSites,cati,strtYrMo,endYrMo)
          #update all outputs here
          dfRes   <- updatedPars$residualData
          dfRes2  <- updatedPars$residualDataNA
          covxx   <- updatedPars$cov
          MVDw    <- updatedPars$mvn
          vpredl  <- updatedPars$pred
          mods    <- updatedPars$newMods 
        }else{message("No univariate outliers")}
      }
    }else{
      missingSites <- setdiff(removeOutlier[[1]], siteOutliers[[1]])
      str2 <- paste(missingSites, collapse= ", ")
      warning(paste0("Outliers were not removed because there exists site(s) ",
                     str2," in removeOutliers that are not in siteOutliers."))}
  }
  if(writeMat){
    write.mat(covxx,filename = "covSites.mat")
  }
  my_list <- list("listMod" = mods, "cov" = covxx, "sites" = cati,
                  "mvn" = MVDw, "univariateTest" = univariateTest, "residualData" = dfRes[,-1],"residualDataNA" = dfRes2[,-1],
                  "rosnerTest" = rosnerT, "pred" = vpredl, "nOutliers"= noutliers)
  return(my_list)
}
