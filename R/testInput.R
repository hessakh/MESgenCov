# Trial input for tests for getCov and getCov2

weeklyB           <- FALSE
startdateStr      <- "01/01/83 00:00"
enddateStr        <- "12/31/90 00:00"
use36             <- TRUE
siteAdd           <- NULL
outliersDates     <- NULL
outlierDatesbySite <- NULL
showOutliers      <- FALSE
siteOutliers      <- NULL
comp              <- c("K","NO3")
plotMulti         <- TRUE
plotB             <- FALSE
sitePlot          <- NULL
plotAll           <- FALSE
writeMat          <- FALSE
seas              <- 12
r                 <- 2
k                 <- 5
p                 <- 1
allComps <- c("Ca","Mg","K","Na","NH4","NO3","Cl","SO4")
perm     <- combn(allComps,2)

cols1 <- dim(perm1)[2]
for (i in perm1){
  comp <- perm[1:2,i]
  print(perm[1:2,i])
  getCov(weeklyB,startdateStr,enddateStr,
        use36,siteAdd,outliersDates,outlierDatesbySite,showOutliers,
        siteOutliers,comp,plotMulti,plotB,sitePlot,plotAll,writeMat,seas,r,k)
}
#1 is decent, 2 is better, 3 is better, 4 worse than 3,5 redo,
#6 is worse than 3, 7 similar to 4, 8 comparable to 3, 9 sim to 4, 10 worth showing prof Lee,
#11 sim to 10,12 flipped S shape, 13-15 sim to 4
#16 bumpy, 17 sim to 4, 18 sim 4?,19 sim 3, 20 bumpy
#21 s-ish close to 3, 22 sim 3 flipped s , 23 decent, 24 just under , 25 sim 23
#26 ok, 27 ok , 28 sim 4
