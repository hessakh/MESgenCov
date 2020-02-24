#' @import lubridate
#' @keywords internal

formatDateTime <- function(data,col){
  data1 <- data
  data1[,col] <- as.POSIXct(data1[,col],
                            format = "%Y-%m-%d %H:%M")
  if(is.na(data1[1,col])){ #try a different format
    data1 <- data
    data1[,col] <- as.POSIXct(data1[,col],
                              format = "%m/%d/%Y %I:%M:%S %p")
    if(is.na(data1[1,col])){
      data1 <- data
      data1[,col] <- as.POSIXct(data1[,col],
                                format = "%m/%d/%y %H:%M")
      if(is.na(data1[1,col])){
        data1 <- data
        data1[,col] <- as.POSIXct(data1[,col],
                                  format = "%Y-%m-%d %H:%M:%S")
        if(is.na(data1[1,col])){#try guess format
          message("still guessing")
          data1 <- data
          #based on 1st date use lubridate to guess column
          gf    <- guess_formats(data1[1,col],c("mdyT","ymdT","dmyT","ymdHM"))
          gf1   <- unname(gf[length(gf)])
          data1[,col] <- as.POSIXct(data1[,col],
                                    format = gf1)
          if(is.na(data1[1,col])){warning("Could not find proper date-time format for data, change source code to include proper format")}
        }
      }
    }
  }
  return(data1)
}
