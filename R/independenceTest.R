#' returns a vector of siteID's of sites with the most data for the specified time period for specific compunds or pH
#' 
#' @param dfRes       Data frame od residuals
#' @return             Returns result of independence test and correlation matrix
#' @export
#' @examples independenceTest(data.frame(matrix(runif(100, min= -1, max =1), ncol=10))
#' @export
independenceTest <- function(dfRes){
  m <- dim(dfRes)[2]
  n <- m*(m-1)/2
  R <- cor(dfRes)
  testStat <- -(n-(2*m+5)/6)*determinant(R, logarithm = TRUE)$modulus
  chisqT  <- qchisq(.95, df=n)
  biChiSq <- testStat <= chisqT
  out <- structure(list(
     u = testStat, chiSq = chisqT, independent = biChiSq)
    ,.Names = c("chisq dist likelihood ratio","chisq","independent") 
    ,row.names = c(NA, -1L) 
    ,class = "data.frame")
    
  result <- list("test" = out, "corr" = R)
}