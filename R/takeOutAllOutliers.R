#'@keywords internal
#' returns LambertW transformed residuals
#' @param dfRes       dataframe, 
#' @param plotMulti   

#' @return             Returns a dataframe of LambertW transformed residuals

takeOutAllOutliers <- function(sites, rosnerResult,cati){

  datesBySiteVec <- NULL
  outSites       <- NULL
  for(i in 1:length(sites)){
    index1 <- match(sites[i],cati)
    stats  <- rosnerResult[[index1]]$all.stats
    outliers <- stats[stats$Outlier == TRUE,]
    if(dim(outliers)[1] > 0){
      outlierst <- outliers$Obs.Num
      outSites  <- c(outSites, sites[i])
      datesBySiteVec <- c(datesBySiteVec, cati[index1])
      for (j in outlierst){
        datesBySiteVec <- c(datesBySiteVec, as.numeric(j))
      }
    }
  }
  list_1 <- list("outBySite" = datesBySiteVec, "sites" = outSites)
  return(list_1)
}