newRes <- NULL
for (i in 1:15){
  # Replicate parts of the analysis in Goerg (2011)
  data("dfRes50")
  y <- dfRes50[,i]
  test_normality(y)
  fit.gmm <- IGMM(y, type = "hh")
  summary(fit.gmm) # gamma is significant and positive
  plot(fit.gmm)
  # Compare empirical to theoretical moments (given parameter estimates)
  # moments.theory <-
  #   mLambertW(theta = list(beta = fit.gmm$tau[c("mu_x", "sigma_x")],
  #                          gamma = fit.gmm$tau["gamma"]),
  #             distname = "normal")
  # TAB <- rbind(unlist(moments.theory),
  #              c(mean(y), sd(y), skewness(y), kurtosis(y)))
  # rownames(TAB) <- c("Theoretical (IGMM)", "Empirical")
  # TAB
  x <- get_input(y, fit.gmm$tau)
  test_normality(x) # input is normal -> fit a Lambert W x Gaussian by MLE
  fit.ml <- MLE_LambertW(y, type = "s", distname = "normal", hessian = TRUE)
  summary(fit.ml)
  plot(fit.ml)
  newRes <- cbind(newRes,x)
}
colnames(newRes) <- cati

#produce univatiate tests for transformed data
MVDw15 <- mvn(newRes, subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "histogram", multivariatePlot = "none", multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", showOutliers = FALSE,showNewData = FALSE)
univariateTest15 <- MVDw15$univariateNormality

#test of data w/o tranformation for comparison
MVDw <- mvn(dfRes50, subset = NULL, mvnTest = "mardia", covariance = TRUE, tol = 1e-25, alpha = 0.5, scale = FALSE, desc = TRUE, transform = "none", univariateTest = "SW",  univariatePlot = "histogram", multivariatePlot = "none", multivariateOutlierMethod = "none", bc = FALSE, bcType = "rounded", showOutliers = FALSE,showNewData = FALSE)
univariateTest <- MVDw$univariateNormality
