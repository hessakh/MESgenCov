#'@keywords internal
combinePrecipitation <- function(mCat,preCSV0,nameFile){
  preCSV <- preCSV0
  #format
  if(!is.null(preCSV)){
    index1 <- match("starttime", colnames(preCSV))
    index2 <- match("endtime"  , colnames(preCSV))
    preCSV <- formatDateTime(preCSV,index1)
    preCSV <- formatDateTime(preCSV,index2)}

  for (i in 1:length(mCat)){
    fileString1 <- paste("NTN-",toString(mCat[i]),"-d.csv", sep = "", collapse = 		NULL)
    tryCatch(
      {tempCSV <- read.csv(fileString1, stringsAsFactors=FALSE)
      on.exit(close)
      index1  <- match("starttime", colnames(tempCSV))
      index2  <- match("endtime"  , colnames(tempCSV))
      tempCSV <- formatDateTime(tempCSV,index1)
      tempCSV <- formatDateTime(tempCSV,index2)
      preCSV  <- rbind(preCSV,tempCSV)
      closeAllConnections() }
      , error = function(e) message(paste("Missing site in dir site",toString(mCat[i])))
    )

  }
  write.csv(preCSV, file = nameFile)
  return(preCSV)
}

