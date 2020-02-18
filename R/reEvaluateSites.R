#' analyzes data without outliers
#' @keywords internal

reEvaluateSites <- function(dfInp,t2,to,startdate,enddate,totT,
                            y0,y,y2,co,inter,e,mods,vpredl, outlierDatesbySite, outSites,cati,strtYrMo,endYrMo,diffYrm){
  df <- dfInp
  weeklyB      <- df$weeklyB
  startdateStr <- df$startdateStr
  enddateStr   <- df$enddateStr
  use36        <- df$use36
  siteAdd      <- df$siteAdd
  outlierDatesbySite <- c(outlierDatesbySite,df$outlierDatesbySite)
  siteOutliers <- df$siteOutliers
  comp         <- df$comp
  plotMulti    <- df$plotMulti
  sitePlot     <- df$sitePlot
  plotAll      <- df$plotAll
  writeMat     <- df$writeMat
  seas         <- df$seas
  r            <- df$r
  k            <- df$k
  
  
  rv    <- integer(5); if(r != 0){rv[1:r] <- rep(1L,r)}
  kv    <- integer(5); if(k != 0){kv[1:k] <- rep(1L,k)}
  kk <- k
  
  obs    <- comp
  si   <- 1
  obsi <- 1
  siObs<- length(obs)
  oc   <- 1

  if(plotAll){
    grDevices::dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
    graphics::par(mfrow = c(2,2),
        oma = c(0,0,0,0) + 0.1,
        mar = c(0,0,0,0) + 0.1)
  }
  for(s in 1:length(outSites)){
    #need to test this part for multiple compounds
    ss <- s # index for site plots, see bottom of loop
    s  <- match(outSites[s],cati)
    si <- s
    s  <- s + (obsi-1)*length(cati)

    sitem <- cbind(t2[[s]], 0, y0[[s]],y2[[s]])
    colnames(sitem) <- c("t","filler","Conc","log")
    sitem <- data.frame(sitem)
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
      t2 <- sitem$t
      sitem     <- stats::na.omit(sitem)
      y1 <- sitem[,4]
      t <- sitem$t
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
      
      #store predicted value vector
      if(length(to[[s]]) < totT){
        new   <- data.frame(1:totT)
        vpred <- predict(mod,newdata = new)
      }else{vpred <- predict(mod)}
      
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
        if(ss%%4 == 1 && ss!=1){
          dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
          graphics::par(mfrow = c(2,2))}
        graphics::par(mar=c(4,4,2,2))
        tc <- 1:totT
        graphics::plot(t,y1, ylab="Log sulfate concentration",main = paste(cati[si],obs[obsi]), xlab = "t (months)")
        graphics::par(new=TRUE)
        graphics::lines(x = tc, y = vpred, col ="blue")
      }
    
    si <- si + 1
  }
  ##################
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

  graphics::par(mar=c(1,1,1,1))
  MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "none", multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", showOutliers = FALSE,showNewData = FALSE)
  univariateTest <- MVDw$univariateNormality
  MVDw
if(plotMulti){
    dev.new(width = 8, height = 5, noRStudioGD = TRUE)
    graphics::par(mfrow=c(1,2))
    MVDw <- mvn(dfRes[,-1], subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", multivariatePlot = "qq", multivariateOutlierMethod = "quan", bc = FALSE, bcType = "rounded", showOutliers = TRUE, showNewData = FALSE)
    MVDw
  }
big_list <- list("mvn" = MVDw, "residualData" = dfRes, "cov" = covxx, "pred" = vpredl)
return(big_list)
}

