## Taken from https://stackoverflow.com/questions/1995933/number-of-months-between-two-dates/1996404

monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01")); 
lt$year*12 + lt$mon } 
