#'creates covariance matrix from normalized NADP monitor data
#'gamma version with optional inputs
#'@import   MVN
#'@import   rmatio
#'@import   EnvStats
#'@param df data frame of input

#'@return   List of model summaries at each site, covariance matrix and plots if inputted as T
#'@export
#'@examples
#'   getCov(structure(list(
#'weeklyB = FALSE, startdateStr = "01/01/83 00:00", enddateStr   = "12/31/86 00:00", 
#'use36 = TRUE, siteAdd = NULL, outlierDatesbySite = NULL, showOutliers = FALSE,
#'siteOutliers = NULL, comp = "SO4", plotMulti = FALSE, plotB = FALSE, sitePlot = NULL,
#'plotAll = FALSE, writeMat = FALSE, seas = 12, r = 1, k = 1)
#',.Names = c("weeklyB","startdateStr","enddateStr","use36","siteAdd",
#'            "outlierDatesbySite","showOutliers","siteOutliers","comp","plotMulti",
#'            "plotB","sitePlot","plotAll","writeMat","seas","r","k") 
#',row.names = c(NA, -1L) 
#',class = "data.frame"))


getCov<- function(df){
  p=1 #for added functionality in the future
  #get data if it's not in the working directory
  if(!exists("weeklyCSV") || !exists("preDailyCSV")){
    try({data("weeklyCSV"); data("preDailyCSV")})
  }
  #check, still doesn't exist?
  if(!exists("weeklyCSV") || !exists("preDailyCSV")){
    stop("Missing files, running code that downloads necessary files from the NADP site")
  }
  #check for multiple pollutants
  if(length(df$comp) > 1){
    stop("Package can only model data from one pollutant/observed variable")
    #later versions of package will be able to handle this
  }
  
  ##code for optional input functionality
  weeklyB      <- df$weeklyB
  startdateStr <- df$startdateStr
  enddateStr   <- df$enddateStr
  use36        <- df$use36
  siteAdd      <- df$siteAdd
  outlierDatesbySite <- df$outlierDatesbySite
  siteOutliers <- df$siteOutliers
  comp         <- df$comp
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
  
  #####
  
  conCSV <- weeklyCSV
  preCSV <- na.omit(preDailyCSV)
  
  #filter out unnecessary columns
  conCSVf  <- conCSV[,-6:-31]
  conCSVf  <- na.omit(conCSVf)
  preCSVf  <- preCSV[,-1]
  
  
  startdate <- as.POSIXct(startdateStr, format = "%m/%d/%y %H:%M")
  enddate   <- as.POSIXct(enddateStr  , format = "%m/%d/%y %H:%M")
  nonneg <- 0
  if(weeklyB){
    totT <- as.integer(difftime(enddate,startdate,units="weeks"))
    nonneg <- 1
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
  
  cati   <- union(cat36,siteAdd[[1]])
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
  preCSVf <- preCSVf[preCSVf$amount>-0.0001,]     #filter out -/ive values
  
  #initialize
  
  rv    <- integer(5); if(r != 0){rv[1:r] <- rep(1L,r)}
  kv    <- integer(5); if(k != 0){kv[1:k] <- rep(1L,k)}
  kk <- k
  
  tl   <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp
  to   <- vector(mode="list", length=(length(cati)*length(obs))) #time stamp
  y0    <- vector(mode="list", length=(length(cati)*length(obs))) #monthly/weekly concentraition values
  y    <- vector(mode="list", length=(length(cati)*length(obs))) #model output for predicted line 	omitted NA values
  y2   <- vector(mode="list", length=(length(cati)*length(obs))) #model output for predicted line w/ NA values
  co   <- vector(mode="list", length=(length(cati)*length(obs))) #model coefficient
  inter<- vector(mode="numeric",length=(length(cati)*length(obs))) #model intercept
  e    <- vector(mode="list", length=(length(cati)*length(obs))) #error
  mods <- vector(mode="list", length=(length(cati)*length(obs))) #model summaries
  vpredl <- vector(mode="list", length=(length(cati)*length(obs))) #list of vector of predictions
  
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
      tl[[s]]     <- NA
      to[[s]]     <- NA
      e[[s]]      <- NA
      y[[s]]      <- NA
      y0[[s]]     <- NA
      y2[[s]]     <- NA
      co[[s]]     <- NA
      inter[s]    <- NA
      mods[[s]]   <- NA
      vpredl[[s]] <- NA
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
      }else{       sitem <- weeklyConc(bi,site,siObs,obs,obsi,startdate)}
      
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
      minsite3 <- min(na.omit(sitem[,3]))
      if(!weeklyB || minsite3 >= 0.0001){
        minsite3 <- 0
        nonneg   <- 0
      }else if(minsite3 == 0){
        nonneg = 0.00001
      }
      sitem[,4] <- log(sitem[,3] + abs(minsite3) + nonneg)
      y0[[s]]   <- sitem[,3]
      y2[[s]]   <- sitem[,4]
      t2        <- sitem$t
      sitem     <- na.omit(sitem)
      colnames(sitem) <- c(cn, "log")
      y1 <- sitem[,4]
      t <- sitem$t
      tsitem <- sitem$t
      cyclicTrend <- (I(cos(t*(2*pi/seas))^p)   + I(sin(t*(2*pi/seas))^p))*kv[1]   +
                     (I(cos(t*(2*pi/seas)*2)^p) + I(sin(t*(2*pi/seas)*2)^p))*kv[2] +
                     (I(cos(t*(2*pi/seas)*3)^p) + I(sin(t*(2*pi/seas)*3)^p))*kv[3] +
                     (I(cos(t*(2*pi/seas)*4)^p) + I(sin(t*(2*pi/seas)*4)^p))*kv[4] +
                     (I(cos(t*(2*pi/seas)*5)^p) + I(sin(t*(2*pi/seas)*5)^p))*kv[5] +
                   t*rv[1] + (t^2)*rv[2] + (t^3)*rv[3] + (t^4)*rv[4] + (t^5)*rv[5]
      df  <- data.frame(cbind(y1,cyclicTrend))
      df  <- data.frame(cbind(df,t))
      mod <- definelm(y1,t,df,r,kk,seas)
      er  <- residuals(mod)
      summary(mod)
      mods[[s]] <- summary(mod)
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
      
      kl    <- 1
      #store model coefficients
      coef <- NULL
      for(j in 1:16){
        if(is.na(unname(mod$coefficients[j]))){
          coef <- c(coef, inter1)
        }else{coef <- c(coef, unname(mod$coefficients[j]))}
      }
      
      #construct predicted value vector
      vpred <- rep(inter1,totT)
      coefi <- 2    # coefficient index
      for (z in 1:totT){
        for (m in 1:kk){
          vpred[z] <- vpred[z] + coef[coefi]*cos(z*(2*pi/seas)*m) + coef[coefi+1]*sin(z*(2*pi/seas)*m)
          coefi <- coefi+2
        }
        coefi = coefi + (5-kk)*2 #jump to time var coef
        for (n in 1:r){
          vpred[z] <- vpred[z] + coef[coefi]*(z^n)
          coefi <- coefi+1
        }
        coefi <- 2 #restart index for new point
      }
      
      for(i in 1:(totT+maxfi-1)){
        if(t[kl] != i-fi){ #if t is skipped
          cy1 <- vpred[i-fi]
          #cy2 + rnorm(1,0,se1) #for maintaining variance
          ry1 <- y1[kl:length(y1)]
          ry2i <- y2i[kl:length(y2i)]
          e0 <- rnorm(1,0,se1)
          y1 <- c(y1[1:kl-1],cy1+e0,ry1)
          y2i <- c(y2i[1:kl-1],NA,ry2i) #missing = 0 rn, rnorm(1,0,se1)
          er <- c(er[1:kl-1],e0,er[kl:length(er)]) # change to stochastic #change to 0
          t  <- c(t[1:kl-1],i-fi,t[kl:length(t)])
          fi <- fi + 1
        }else{kl = kl + 1}
      }
      if(length(t)>totT){t<-t[1:totT];y2i<-y2i[1:totT]; y1 <- y1[1:totT];er<-er[1:totT]}
      
      #store parameters
      tl[[s]]    <- t
      e[[s]]     <- unname(er)
      y[[s]]     <- y1
      #y2[[s]]    <- y2i
      co[[s]]    <- coef
      inter[s]   <- inter1
      vpredl[[s]] <- vpred
      #plot all sites if indicated by user
      if(plotAll == T){
        if(s%%4 == 1 && s!=1){
          dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
          par(mfrow = c(2,2))}
        par(mar=c(4,4,2,2))
        tc <- 1:totT
        plot(t2,y2[[s]],ylim=c(-2, 3), ylab="Log sulfate concentration",main = paste(cati[si],obs[obsi]), xlab = "t (months)")
        par(new=TRUE)
        lines(x = tc, y = vpred, col ="blue")
        #plot(t, vpred, type="l",col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
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
  par(mar=c(1,1,1,1))
  MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "none", multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", showOutliers = FALSE,showNewData = FALSE)
  univariateTest <- MVDw$univariateNormality
  MVDw
  if(plotMulti){
    dev.new(width = 8, height = 5, noRStudioGD = TRUE)
    par(mfrow=c(1,2))
    MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "qq", multivariateOutlierMethod = "quan", bc = FALSE, bcType = "rounded", showOutliers = TRUE, showNewData = FALSE)
    MVDw
  }
  siteRosner = NULL
  if(showOutliers == TRUE){
    #find outliers in each site
    #Store Site outliers
    i <- match(siteOutliers,colnames(dfRes))
    if(is.na(i)){
      message("siteOutliers string doesn't match any colnames, these are the options:")
      print(colnames(dfRes))
      siteRosner = NULL}
    else{siteRosner = rosnerTest(dfRes[,i])}
    
  }
  if(plotB == T){
    dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
    par(mar=c(4,4,2,2))
    i = match(sitePlot[1], cati)
    tc <- 1:totT
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration", xlab = "t (months)", main = toString(sitePlot[1]))
    par(new=TRUE)
    plot(tc, vpredl[[i]], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "",xlab="")
  }
  if(writeMat){
    write.mat(covxx,filename = "covSites.mat")
  }
  my_list <- list("listMod" = mods, "cov" = covxx, "sites" = cati,
                  "mvn" = MVDw, "univariateTest" = univariateTest, "residualData" = dfRes[,-1],
                  "rosnerTest" = siteRosner, "pred" = vpredl)
  return(my_list)
}
