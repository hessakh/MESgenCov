#'@keywords internal
## Code by Dirk Eddelbuettel taken from 
##https://stackoverflow.com/questions/1995933/number-of-months-between-two-dates/1996404

# compute a month difference as a difference between two monnb's
mondf <- function(d1, d2) {monnb(d2) - monnb(d1)}