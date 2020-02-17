#' returns LambertW transformed residuals
#' @keywords internal
#' @importFrom MVN mvn
#' @import LambertW
#' @param dfRes       dataframe, 

#' @return             Returns a dataframe of LambertW transformed residuals
#' @export
#' @examples lambertWfit(data.frame(matrix(runif(n = 100), ncol = 5 )))

lambertWfit <- function (dfRes){
  
  storeType <- NULL
  xhh       <- NULL
  fitgmmhh  <- NULL
  xh        <- NULL
  fitgmmh   <- NULL
  xs        <- NULL
  fitgmms   <- NULL
  newRes    <- dfRes[,1]
  concdf    <- dfRes
  nosites   <- dim(dfRes)[2]-1
  totT      <- dim(dfRes)[1]
  options(warn=-1)
  for (i in 1:nosites){
    types1    <- NULL
    # Replicate parts of the analysis in Goerg (2011)
    yL <- na.omit(concdf[,i+1])
    try(fitgmmhh <- IGMM(yL, type = "hh"),silent = TRUE)
    if(!is.null(fitgmmhh)){
      xhh <- get_input(yL, fitgmmhh$tau)
      temphh <- shapiro.test(xhh)
      pxhh   <- temphh$p.value 
      names(pxhh) <- "hh"
      types1 <- c(types1, "hh")
    }else{pxhh <- NULL}
    try(fitgmmh <- IGMM(yL, type = "h"),silent = TRUE)
    if(!is.null(fitgmmh)){
      xh <- get_input(yL, fitgmmh$tau)
      temph  <- shapiro.test(xhh)
      pxh    <- temph$p.value 
      names(pxh) <- "h"
      types1 <- c(types1, "h")
    }else{pxh <- NULL}
    try(fitgmms <- IGMM(yL, type = "s"),silent = TRUE)
    if(!is.null(fitgmms)){
      xs <- get_input(yL, fitgmms$tau)
      temps  <- shapiro.test(xs)
      pxs <- temps$p.value
      names(pxs) <- "s"
      types1 <- c(types1, "s")
    }else{pxs <- NULL}
    
    #x <- cbind(xhh,xh,xs)
    pvals <- c(pxhh,pxh,pxs)
    maxT <- types1[which.max(pvals)]
    #x <- x[,which.max(pvals)]
    storeType <- c(storeType,maxT)
    
    ###### get input after filling in missing values
    if(length(yL) < totT){
      if(maxT == "hh"){
        fitgmm <- fitgmmhh
      }else if(maxT == "h"){
        fitgmm <- fitgmmh
        #print(paste0("i = ", i))
      }else{
        fitgmm <- fitgmms
      }
      newResi <- fillInLambert(concdf[,i+1],fitgmm, maxT) 
    }else{newResi <- yL}
    #print(mean(x))
    newRes <-  cbind(newRes,newResi)
    fitgmm <- NULL
    xhh    <- NULL
    xh     <- NULL
    xs     <- NULL
  }
  options(warn=0)
  newRes <- as.data.frame(newRes)
  colnames(newRes) <- colnames(dfRes)
  
  my_list <-  newRes
  return(my_list)
}


