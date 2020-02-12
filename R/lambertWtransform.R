
#' returns LambertW transformed residuals
#' @import LambertW
#' @import MVN
#' @param dfRes       dataframe, 
#' @param plotMulti   Binary, 

#' @return             Returns a dataframe of LambertW transformed residuals
#' @export
#' @examples lambertWtransform(data.frame(matrix(runif(n = 100), ncol = 5 )),plotMulti = FALSE, writeMat = FALSE)

lambertWtransform <- function (dfRes,plotMulti,writeMat){

  storeType <- NULL
  xhh       <- NULL
  fitgmmhh  <- NULL
  xh        <- NULL
  fitgmmh   <- NULL
  xs        <- NULL
  fitgmms   <- NULL
  newRes    <- NULL
  concdf    <- dfRes
  nosites   <- dim(dfRes)[2]
  totT      <- dim(dfRes)[1]
  types1    <- c("hh","h","s")
  options(warn=-1)
  for (i in 1:nosites){
    # Replicate parts of the analysis in Goerg (2011)
    yL <- na.omit(concdf[,i])
    try(fitgmmhh <- IGMM(yL, type = "hh"),silent = TRUE)
    if(!is.null(fitgmmhh)){
      xhh <- get_input(yL, fitgmmhh$tau)
      temphh <- shapiro.test(xhh)
      pxhh   <- temphh$p.value 
      names(pxhh) <- "hh"
    }else{pxhh <- NULL}
    try(fitgmmh <- IGMM(yL, type = "h"),silent = TRUE)
    if(!is.null(fitgmmh)){
      xh <- get_input(yL, fitgmmh$tau)
      temph  <- shapiro.test(xhh)
      pxh    <- temph$p.value 
      names(pxh) <- "h"
    }else{pxh <- NULL}
    try(fitgmms <- IGMM(yL, type = "s"),silent = TRUE)
    if(!is.null(fitgmms)){
      xs <- get_input(yL, fitgmms$tau)
      temps  <- shapiro.test(xs)
      pxs <- temps$p.value
      names(pxs) <- "s"
    }else{pxs <- NULL}

    x <- cbind(xhh,xh,xs)
    pvals <- c(pxhh,pxh,pxs)
    maxT <- types1[which.max(pvals)]
    x <- x[,which.max(pvals)]
    storeType <- c(storeType,maxT)

    ###### get input after filling in missing values
    if(length(yL) < totT){
      if(maxT == "hh"){
        fitgmm <- fitgmmhh
      }else if(maxT == "h"){
        fitgmm <- fitgmmh
      }else{
        fitgmm <- fitgmms
      }
      newResi <- fillInLambert(concdf[,i],fitgmm, maxT) 
      x       <- get_input(newResi, fitgmm$tau)
    }
    #print(mean(x))
    x <- x - mean(x)
    newRes <-  cbind(newRes,x)
    fitgmm <- NULL
    xhh    <- NULL
    xh     <- NULL
    xs     <- NULL
    # if(i %% 4 == 1){
    #   dev.new()
    #   par(mfrow = c(2,2),
    #       oma = c(0,0,0,0) + 0.1,
    #       mar = c(4,4,0.5,0.5) + 0.1)
    # }
    # hist(newResi)
    # hist(concdf[,i])
  }
  options(warn=0)
  newRes <- as.data.frame(newRes)
  colnames(newRes) <- colnames(dfRes)
  
  #produce univatiate tests for transformed data
  if(plotMulti){
    plotM = "qq"
    dev.new(width = 8, height = 5, noRStudioGD = TRUE)
    par(mfrow = c(1,2))
  }else{plotM = "none"}
  MVDw15 <- mvn(newRes, subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, 
                scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "none", 
                multivariatePlot = plotM, multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", 
                showOutliers = FALSE,showNewData = FALSE)
  covxx <- data.frame(cov(newRes))
  #corrxx <- data.frame(cor(newRes))
  write.mat(corrxx,filename = "corrSitesLambertW1.mat")
  if(writeMat){
    write.mat(covxx,filename = "covSitesLambertW.mat")
  }

  
  my_list <- list("newResiduals" = newRes, "mvn" = MVDw15, "univariateTest" =  MVDw15$univariateNormality,
                  "cov" = covxx)
  return(my_list)
}


