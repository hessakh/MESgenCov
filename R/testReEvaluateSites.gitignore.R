#for testing reEvaluate sites

# maxds5  <- maxDistSites("01/01/86 00:00","12/31/88 00:00",5,50,"SO4",1) 
# 
# ##create input dataframe
# df50 <- defaultInput
# rand <- maxds5
# df50$siteAdd <- list(rand$finalList)
# df50$startdateStr <- rand$startDate
# df50$use36 <- FALSE
# df50$weeklyB <- FALSE
# df50$comp <-  rand$comp
# df50$plotAll <-  T
# df50$enddateStr <- rand$endDate
# df50$writeMat   <- F
# df50$k          <- 3
# df50$r          <- 1
# df50$siteOutliers   <- list(rand$finalList)
# df50$removeOutliers <- list(rand$finalList)
# test53 <- getCov(df50)
# test53$mvn$multivariateNormality
