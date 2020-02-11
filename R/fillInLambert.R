#'@keywords internal
#'This takes a fitted Lambert distribution, residuals and a type of LambertW as input
#' and outputs new residuals with missing values sampled from the LambertW fitted distribution
# input: fitgmm, residual data
#output: filled in residuals

fillInLambert <- function(resNA,fitgmm,maxT){
  newRes <- NULL
  noNA <- length(resNA) - length(na.omit(resNA))
  if(maxT == "hh"){
    g <- fitgmm$skewness.x
    t <- fitgmm$tau.init
    d <- c(fitgmm$tau.init[3:4])
    b <- c(fitgmm$tau.init[1:2])
    set.seed(i+100)
    estimate <- rLambertW(n=noNA, distname = "normal", theta =  list(beta = b, gamma = 0, delta = d), tau = t)
  }else if(maxT == "h"){
    g <- fitgmm$skewness.x
    t <- fitgmm$tau.init
    d <- c(fitgmm$tau.init[3:4])
    b <- c(fitgmm$tau.init[1:2])
    set.seed(i+100)
    estimate <- rLambertW(n = noNA, distname = "normal", theta =  list(beta = b, gamma = g), delta = d, tau = t)
  }else{
    g <- fitgmm$tau.init[3]
    t <- fitgmm$tau.init
    b <- c(fitgmm$tau.init[1:2])
    set.seed(i)
    estimate <- rLambertW(n = noNA, distname = "normal", theta = list(beta = b, gamma = g))
  }
  j <- 1
  for (i in 1:length(resNA)){
    if(is.na(resNA[i])){
      resNA[i] <- estimate[j]
      j <- j + 1
    }
  }
  return(resNA)
}
