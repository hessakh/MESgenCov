plotOpt <- function(n,sites,y,t){ #not sure how to include comps
  if(length(sitePlot)==4){
    dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
    par(mfrow=c(2,2))
    #plot logConc vs t and model
    par(mar=c(4,4,2,2))
    i = match(sitePlot[1], cati)  #add string siteID here for desired site logConc plot, change plot title accordingly
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration", xlab = "t (months)", main = cati[i])
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(mar=c(4,4,2,2))
    i = match(sitePlot[2], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration",main = cati[i], xlab = "t (months)")
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(new=F)
    par(mar=c(4,4,2,2))
    i = match(sitePlot[3], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration",xlab = "t (months)", main = cati[i])
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(mar=c(4,4,2,2))
    i = match(sitePlot[4], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration",xlab = "t (months)", main = cati[i])
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
  }else if(length(sitePlot)==3){
    dev.new(width = 6, height = 5.5, noRStudioGD = T, unit = "in")
    par(mfrow=c(2,2))
    #plot logConc vs t and model
    par(mar=c(4,4,2,2))
    i = match(sitePlot[1], cati)  #add string siteID here for desired site logConc plot, change plot title accordingly
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration", xlab = "t (months)", main = cati[i])
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(mar=c(4,4,2,2))
    i = match(sitePlot[2], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration",main = cati[i], xlab = "t (months)")
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(new=F)
    par(mar=c(4,4,2,2))
    i = match(sitePlot[3], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration",xlab = "t (months)", main = cati[i])
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
  }else if(length(sitePlot)==2){
    dev.new(width = 8, height = 5, noRStudioGD = T, unit = "in")
    par(mfrow=c(2,2))
    #plot logConc vs t and model
    par(mar=c(4,4,2,2))
    i = match(sitePlot[1], cati)  #add string siteID here for desired site logConc plot, change plot title accordingly
    print(i)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration", xlab = "t (months)", main = cati[i])
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(mar=c(4,4,2,2))
    i = match(sitePlot[2], cati)
    plot(tc,y2[[i]],ylim=c(-2, 3), ylab="Log sulfate concentration",main = cati[i], xlab = "t (months)")
    par(new=TRUE)
    plot(tc,co[[i]][1]*I(cos(tc*(2*pi/seas))^p)*kv[1]   + co[[i]][2]*I(cos(tc*(2*pi/seas)*2)^p)*kv[2] +
           co[[i]][3]*I(cos(tc*(2*pi/seas)*3)^p)*kv[3] + co[[i]][4]*I(cos(tc*(2*pi/seas)*4)^p)*kv[4] +
           co[[i]][5]*I(cos(tc*(2*pi/seas)*5)^p)*kv[5] + co[[i]][6]*tc*rv[1] + co[[i]][7]*(tc^2)*rv[2] +co[[i]]					[8]*(tc^3)*rv[3]+
           co[[i]][9]*(tc^4)*rv[4] + co[[i]][10]*(tc^5)*rv[5] + inter[i], type="l"
         ,col ="blue",ylim=c(-2, 3),ylab = "", xlab ="")
    par(new=F)
    par(mar=c(4,4,2,2))
  }
}
