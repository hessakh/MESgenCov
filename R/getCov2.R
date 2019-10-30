#'creates covariance matrix from normalized NADP monitor data
#'getCov2 uses the old model
#'@import MVN
#'@import rmatio
#'@import EnvStats
#'@param weeklyB        Binary, analyze weekly data, if F monthly data is analyzed
#'@param startdateStr   String, when to start analyzing data, format = "m/d/y H:M"
#'@param enddateStr     String, when to stop analyzing data, format = "m/d/y H:M"
#'@param use36          Binary, use 36 default sites
#'@param siteAdd         Vector, sites to add
#'@param outliersDates  Integer, list of months/weeks (integers, 1 => 1st month/week) to omit
#'@param outlierDatesbySite Vector of String, Integer, list of sites and months to omit
#'@param showOutliers   Binary, show outliers for specific site
#'@param siteOutliers   String, specify siteID to check outliers
#'@param comp           Vector of Strings, compunds to include in analysis
#'@param plotMulti      Binary, plot multivariate analysis qq plot and outlier test
#'@param plotB          Binary, display plots of fit
#'@param sitePlot       Vector of Strings, of sites to plot
#'@param plotAll        Binary, plot all sites
#'@param writeMat       Binary, download cov.mat file in working directory
#'@param seas           Integer, period for model
#'@param r              Integer, degree of polynomial to fit
#'@param k              Integer, stopping term for fourier series
#'@param p              Integer, power on trignometric terms
#'@return               List of model summaries at each site, covariance matrix and plots if inputted as T
#'@export
#'@examples
#'   getCov2(FALSE,"01/01/83 00:00","12/31/86 00:00",TRUE,
#'   c("AK01","PA02"),c(5,10),c("IL19",9, "VT01",1,33),TRUE,
#'   "VT01",c("SO4"),FALSE,FALSE,NULL,FALSE,FALSE,12,1,3,1)


getCov2 <- function(weeklyB,startdateStr,enddateStr,
                   use36,siteAdd,outliersDates,outlierDatesbySite,showOutliers,
                   siteOutliers,comp,plotMulti,plotB,sitePlot,plotAll,writeMat,seas,r,k,p){

  #get data if it's not in the working directory
  if(!exists("weeklyCSV") || !exists("preDailyCSV")){
    try({data("weeklyCSV"); data("preDailyCSV")})
  }
  #check, still doesn't exist?
  if(!exists("weeklyCSV") || !exists("preDailyCSV")){
    message("Missing files, running code that downloads necessary files from the NADP site")
    getDataOffSite()
  }

  conCSV <- weeklyCSV
  preCSV <- na.omit(preDailyCSV)

  #filter out unnecessary columns
  conCSVf  <- conCSV[,-6:-31]
  conCSVf  <- na.omit(conCSVf)
  preCSVf  <- preCSV[,-1]


  startdate <- as.POSIXct(startdateStr, format = "%m/%d/%y %H:%M")
  enddate   <- as.POSIXct(enddateStr  , format = "%m/%d/%y %H:%M")
  if(weeklyB){
    totT <- as.integer(difftime(enddate,startdate,units="weeks"))
    if(is.null(seas)){seas = 52}
  }else if(!(floor(mondf(startdate,enddate))==(mondf(startdate,enddate)))){# won't happen because code wont produce decimals
    stop("Number of months is not integer")
  }else{
    totT <- (mondf(startdate,enddate))+1
    if(is.null(seas)){seas = 12}
  }

  if(use36 == T){
    cat36 <- c("AL10","IL11","IL18","IL19","IL35","IL47","IL63","IN34",
               "IN41","MA01","MA13","MD13","MI09","MI26","MI53",
               "NC03","NC34","NC41","NJ99","NY08","NY10","NY20",
               "NY52","NY65","OH17","OH49","OH71","PA15","PA29",
               "PA42","TN00","TN11","VA13","VT01","WI28","WV18")#from Guttorp, Le 1992
  }else{cat36 <- NULL}

  cati   <- union(cat36,siteAdd)
  obs    <- comp

  #add back desired columns based on input in obs e.g add back SO4, pH, NO3 etc.
  for (i in 1:length(obs)){
    obsiCSV           <- match(obs[i],colnames(conCSV))
    cn                <- colnames(conCSVf)
    conCSVf[,5+i]     <- conCSV[,obsiCSV] #
    colnames(conCSVf) <- c(cn, obs[i])
  }

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


  preCSVf$amount    <- as.numeric(preCSVf$amount) #change data type of column
  preCSVf <- preCSVf[preCSVf$amount>-0.0001,]   #filter out -/ive values

  #initialize

  rv    <- integer(5); if(r != 0){rv[1:r] <- rep(1L,r)}
  kv    <- integer(5); if(k != 0){kv[1:k] <- rep(1L,k)}

  tl   <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp
  to   <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp
  y    <- vector(mode="list", length=(length(cati)*length(obs))) #model output for predicted line 	omitted NA values
  y2   <- vector(mode="list", length=(length(cati)*length(obs))) #model output for predicted line w/ NA values
  co   <- vector(mode="list", length=(length(cati)*length(obs))) #model coefficient
  inter<- vector(mode="numeric", length=(length(cati)*length(obs))) #model intercept
  e    <- vector(mode="list", length=(length(cati)*length(obs))) #error
  mods <- vector(mode="list", length=(length(cati)*length(obs))) #error

  betaCos <- NULL
  betaSin <- NULL
  betaT   <- NULL
  inter 	<- NULL
  tl   <- vector(mode="list", length=length(cat)) #time stamp
  y    <- vector(mode="list", length=length(cat)) #model output for predicted line
  e    <- vector(mode="list", length=length(cat)) #error
  mode <- vector(mode="list", length=length(cat)) #stores model

  si   <- 1
  obsi <- 1
  siObs<- length(obs)
  oc   <- 1


  #par(mfrow=c(2,2))
  par(mfrow = c(2,2),
      oma = c(0,0,0,0) + 0.1,
      mar = c(0,0,0,0) + 0.1)
  for(s in 1:(length(cati)*length(obs))){
    #change obs
    if (s%%(length(cati)) == 1 && s != 1){obsi = obsi + 1; si = 1}

    #filter site (weekly concentration data)

    conCSVk <- conCSVf[conCSVf$siteID == toString(cati[si]),]
    preCSVk <- preCSVf[preCSVf$siteID == toString(cati[si]),]

    if(nrow(conCSVk)==0 || nrow(preCSVk)==0){
      #store parameters
      tl[[s]] <- NA
      to[[s]] <- NA
      e[[s]]  <- NA
      y[[s]]  <- NA
      y2[[s]] <- NA
      co[[s]] <- NA
      inter[s]<- NA
      mods[[s]] <- NA
      message(paste0("Missing data for ", cati[si]," ",obs[obsi],
                     " check if site has data for inputted dates in data file weeklyCSV, preDailyCSV"))
    }else{
      #filter out missing data, -/ive values, and order dataframes by date
      site <- conCSVk
      site <- site[site[,5 + obsi]>-0.0001,]
      site <- site[order(site$dateon),]
      preCSVk <- preCSVk[order(preCSVk$starttime,decreasing=F),]
      conCSV$siteID <- toupper(conCSV$siteID) #combine data
      preCSV$siteID <- toupper(preCSV$siteID)

      #aaggregates weekly precipitation
      if (obs[obsi] == "ph"){bi = 1}else{bi = 0}
      site <- appendPre(site,preCSVk,obs,obsi,siObs,bi)
      site <- na.omit(site)

      #aggregate precipitation data monthly and get concentration values
      if(!weeklyB){sitem <- aggregateMonthly(bi,site,siObs,obs,obsi,totT,strtYrMo)
      }else{     sitem <- weeklyConc(bi,site,siObs,obs,obsi,startdate)}
      sitem <- sitem[sitem[,3]>= 0.0000001,]
      #take out outliers
      if (length(outliersDates) != 0){
        outliersDates <- sort(outliersDates, decreasing = F)
        for(i in 1:length(outliersDates)){
          j     <- outliersDates[i]
          index <- match(j, sitem$t)
          if (!is.na(index)){
            sitem[index,3] <-  NA
          }
        }#endfor
      }#endif
      #still need to debug this part
      if (length(outlierDatesbySite) != 0){
        if(outlierDatesbySite[oc] == cati[si]){
          for(i in (oc+1):length(outlierDatesbySite)){
            ifInt  <- suppressWarnings(as.integer(outlierDatesbySite[i]))
            if (!is.na(ifInt)){
              j     <- as.integer(outlierDatesbySite[i])
              index <- match(j, sitem$t)
              if (!is.na(index)){
                sitem[index,3] <- NA
              }
            }else{break}
          }
          oc = i #outlier site counter #skips t
        }
      }

      #deterministic trend for one site
      cn <- colnames(sitem)
      sitem[,4] <- log(sitem[,3])
      colnames(sitem) <- c(cn, "log")
      y2[[s]] <- sitem[,4]
      sitem <- na.omit(sitem)
      y1 <- sitem[,4]
      t <- sitem$t
      cyclicTrend <- cos(t*(pi/6)) + sin(t*(pi/6)) + t
      df <- data.frame(cbind(y1,cyclicTrend))
      df <- data.frame(cbind(df,t))
      mod <- lm(y1 ~ cos(t*(pi/6)) + sin(t*(pi/6)) + t, data = df)
      er  <- residuals(mod)
      summary(mod)
      to[[s]] <- t

      #fill in missing t

      inter1 <- unname(mod$coefficients[1])
      se1  <- sqrt(deviance(mod)/df.residual(mod))
      fi   <- 0 #to loop back
      maxfi<- totT - length(t)
      fe   <- 1:(maxfi) #fake end to keep loop going #ignore if t=48
      t    <- c(t,fe)
      y1   <- c(y1,fe)
      y2i  <- c(y1,fe)
      er   <- c(er,fe)
      mods[[s]] <- summary(mod)
      kl    <- 1

      coef <- NULL
      for(j in 2:11){
        if(is.na(unname(mod$coefficients[j]))){
          coef <- c(coef, 0)
        }else{coef <- c(coef, unname(mod$coefficients[j]))}
      }
      betaCos1 <- unname(mod$coefficients[2]) #store parameters
      betaSin1 <- unname(mod$coefficients[3])
      betaT1   <- unname(mod$coefficients[4])
      inter1 	<- unname(mod$coefficients[1])
      se      <-sqrt(deviance(mod)/df.residual(mod))
      fi <- 0 #to loop back
      maxfi <- totT - length(t)
      fe <- 1:(maxfi) #fake end to keep loop going
      t <- c(t,fe)
      y1 <- c(y1,fe)
      er <- c(er,fe)
      mode <- c(mode,mod)
      k <- 1

      for(i in 1:(totT+maxfi-1)){
        if(t[k] != i-fi){ #if t is skipped
          ry1 <- y1[k:length(y1)]
          cy1 <- inter1 + betaCos1*cos((i-fi)*(pi/6)) + betaSin1*sin((i-fi)*(pi/6)) + (i-fi)*betaT1 + rnorm(1,0,se)
          cy2 <- inter1 + betaCos1*cos((i-fi)*(pi/6)) + betaSin1*sin((i-fi)*(pi/6)) + (i-fi)*betaT1
          y1 <- c(y1[1:k-1],cy1,ry1)
          er <- c(er[1:k-1],cy1-cy2,er[k:length(er)]) #change back to -> cy1-cy2
          t  <- c(t[1:k-1],i-fi,t[k:length(t)])
          fi <- fi + 1
        }else{k = k + 1}
      }
      if(length(t)>totT){t<-t[1:totT];y1 <- y1[1:totT];er<-er[1:totT]}

      #store parameters old mod
      betaCos <- c(betaCos,betaCos1)
      betaSin <- c(betaSin ,betaSin1)
      betaT   <- c(betaT ,betaT1)
      inter 	<- c(inter ,inter1)
      tl[[s]] <- t
      e[[s]]  <- er
      y[[s]]  <- y1
      tc = 1:totT
      if(plotAll == T){
        if(s%%4 == 1 && s!=1){
          dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
          par(mfrow = c(2,2))}
        par(mar=c(4,4,2,2))
        plot(tc,y2[[s]],ylim=c(-2, 3), ylab="Log sulfate concentration",main = paste(cati[si],obs[obsi]), xlab = "t (months)")
        par(new=TRUE)
        plot(tc,  betaCos1*cos(tc*(pi/6)) + betaSin1*sin(tc*(pi/6)) + betaT1*tc + inter1, type="l",col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
      }
    }
    si <- si + 1
  }
  ############################END BIG LOOP##############################

  options(warn=-1)
  dfRes <- NULL
  si    <- 2
  obsi  <- 1
  tc <- 1:totT
  #loop to create dataframe of all residuals
  dfRes <- cbind(tc,e[[1]])
  colnames(dfRes) <- c("t",paste(cati[1],obs[1], sep=""))
  for( i in 2:(length(cati)*length(obs))){
    if (i%%(length(cati)) == 1){obsi = obsi + 1; si = 1}
    if(!is.na(e[[i]])){
      currRes <- cbind(tc,e[[i]])
      colnames(currRes) <- c("t",paste0(cati[si],obs[obsi]))
      dfRes   <- merge(dfRes,currRes, by = "t", all = F)
    }
    si <- si+1
  }
  #create covariance matrix
  covxx <- data.frame(cov(dfRes[,-1]))
  options(warn=0)


  # filename <- paste(nameNewCov, ".csv", sep = "")
  # cname = cname+1
  # write.csv(covxx, file = filename)
  MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "none", multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", showOutliers = FALSE,showNewData = FALSE)
  univariateTest <- MVDw$univariateNormality
  if(plotMulti){
    dev.new(width = 8, height = 5, noRStudioGD = TRUE)
    par(mfrow=c(1,2))
    MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "qq", multivariateOutlierMethod = "quan", bc = FALSE, bcType = "rounded", showOutliers = TRUE, showNewData = FALSE)
    MVDw
  }

  if(showOutliers == TRUE){
    #find outliers in each site
    i = match(siteOutliers, colnames(dfRes))
    rosnerTest(dfRes[,i])
  }
  if(plotB == T){
    dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
    par(mar=c(4,4,2,2))
    i = match(sitePlot[1], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration", xlab = "t (months)", main = toString(sitePlot[1]))
    par(new=TRUE)
    plot(tc,betaCos[i]*cos(tc*(pi/6)) + betaSin[i]*sin(tc*(pi/6)) + betaT[i]*tc + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "",xlab="")
  }
  if(writeMat){
    write.mat(covxx,filename = "covSites.mat")
  }
  my_list <- list("listMod" = mods, "cov" = covxx, "sites" = cati, "mvn"=MVDw, "univariateTest"=univariateTest)
  return(my_list)
}
