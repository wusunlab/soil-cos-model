fileName <- './input/soil_flux_regular_meas.csv'
okFlux <- read.csv( fileName, skip = 0, header = TRUE)

expFunc <- function(T, p) {
  p[1] + p[2] * exp( p[3] * T ) 
}

ArrFunc <- function(T_k, p) {
  p[1] + p[2] * exp( -p[3] / 8.3145 * ((1/T_k) - (1/298.15))) 
}

resd <- function(T,obs,p) {
  obs - p[1] + p[2] * exp( p[3] * T )
}

fit_flux_temp <- list()
fit_flux_temp_r2 <- c()
fit_flux_temp_se <- c()
fit_flux_temp_mae <- c()
Q10 <- c()
Q10_lower <- c()
Q10_upper <- c()

# create an index for soil moisture bins
okFlux$swc_bin <- NaN
for ( i in floor(min(okFlux$soil_moist_Avg, na.rm=TRUE)/0.05):
        floor(max(okFlux$soil_moist_Avg, na.rm=TRUE)/0.05) ) {
  okFlux$swc_bin[okFlux$soil_moist_Avg >= 0.05*i & okFlux$soil_moist_Avg < 0.05*(i+1)] <- i
#   x <- okFlux$soil_temp_C_Avg[okFlux$swc_bin==i]
#   y <- okFlux$fCOS_corr[okFlux$swc_bin==i]
#   sym_size <- okFlux$soil_moist_Avg[okFlux$swc_bin==i]
#   plot(x,y)
#   symbols(x,y,circles=sym_size, inches=1/16)
#   fit_flux_temp[[i]] <- nls(y ~ expFunc(x, p), alg='port',
#                             start=list( p = c(-2, 4, 0.05) ), trace=TRUE,
#                             control=list(tol=1e-6, maxiter = 1000))
}

# create an index for growing stages: growing 1, senescence 2, post-harvest 3
okFlux$stage <- NaN
okFlux$stage[okFlux$doy_local < 120.] <- 1
okFlux$stage[okFlux$doy_local > 120. & okFlux$doy_local < 145. ] <- 2
okFlux$stage[okFlux$doy_local > 145.] <- 3

symbols(okFlux$soil_temp_C_Avg[okFlux$stage==1], okFlux$fCOS_corr[okFlux$stage==1], circles=okFlux$soil_moist_Avg[okFlux$stage==1], inches=1/16)
symbols(okFlux$soil_temp_C_Avg[okFlux$stage==2], okFlux$fCOS_corr[okFlux$stage==2], circles=okFlux$soil_moist_Avg[okFlux$stage==2], inches=1/16)
symbols(okFlux$soil_temp_C_Avg[okFlux$stage==3], okFlux$fCOS_corr[okFlux$stage==3], circles=okFlux$soil_moist_Avg[okFlux$stage==3], inches=1/16)

# create a fitting index: 1 - growing season, 2 - senescence, high SWC, 3 - senescence, low SWC
# 4 - post-harvest
okFlux$fit_ind <- NaN
okFlux$fit_ind[okFlux$doy_local < 120.] <- 1
okFlux$fit_ind[okFlux$doy_local > 120. & okFlux$doy_local < 145. & okFlux$soil_moist_Avg >= 0.2] <- 2
okFlux$fit_ind[okFlux$doy_local > 120. & okFlux$doy_local < 145. & okFlux$soil_moist_Avg < 0.2] <- 3
okFlux$fit_ind[okFlux$doy_local > 145.] <- 4

for (i in 1:4) {
  x <- okFlux$soil_temp_C_Avg[okFlux$fit_ind==i]
  y <- okFlux$fCOS_corr[okFlux$fit_ind==i]
  fit_flux_temp[[i]] <- nls(y ~ expFunc(x, p), alg='default',
                            start=list( p = c(-2, 4, 0.05) ), trace=TRUE,
                            control=list(tol=1e-6, maxiter = 1000))
  fit_flux_temp_r2[i] <- cor(y[is.finite(x) & is.finite(y)],predict(fit_flux_temp[[i]]))^2
  fit_flux_temp_se[i] <- summary(fit_flux_temp[[i]])$sigma
  fit_flux_temp_mae[i] <- sum(abs(residuals(fit_flux_temp[[i]])))/summary(fit_flux_temp[[i]])$df[2]
  Q10[i] <- exp(coef(fit_flux_temp[[i]])[3]*10)
  Q10_lower[i] <- exp((coef(summary(fit_flux_temp[[i]]))[3,1] - coef(summary(fit_flux_temp[[i]]))[3,2])*10)
  Q10_upper[i] <- exp((coef(summary(fit_flux_temp[[i]]))[3,1] + coef(summary(fit_flux_temp[[i]]))[3,2])*10)
}

# activation energy fitting
fit_flux_Ea <- list()
fit_flux_Ea_r2 <- c()
fit_flux_Ea_se <- c()
fit_flux_Ea_mae <- c()
E_act <- c()
E_act_sd <- c()

for (i in 1:4) {
  x <- okFlux$soil_temp_C_Avg[okFlux$fit_ind==i] + 273.15
  y <- okFlux$fCOS_corr[okFlux$fit_ind==i]
  fit_flux_Ea[[i]] <- nls(y ~ ArrFunc(x, p), alg='port',
                            start=list( p = c(0, 2, 30000) ), trace=TRUE,
                            control=list(tol=1e-6, maxiter = 1000))
  fit_flux_Ea_r2[i] <- cor(y[is.finite(x) & is.finite(y)],predict(fit_flux_Ea[[i]]))^2
  fit_flux_Ea_se[i] <- summary(fit_flux_Ea[[i]])$sigma
  fit_flux_Ea_mae[i] <- sum(abs(residuals(fit_flux_Ea[[i]])))/summary(fit_flux_Ea[[i]])$df[2]
  E_act[i] <- coef(fit_flux_Ea[[i]])[3]
  E_act_sd[i] <- coef(summary(fit_flux_Ea[[i]]))[3,2]
}

# raw plot of f_COS vs soil temp
plot(okFlux$soil_temp_C_Avg, okFlux$fCOS_corr)
plot(okFlux$soil_temp_C_Avg, okFlux$fCO2_flow)

# plot exponential fittings
plot(okFlux$soil_temp_C_Avg[okFlux$fit_ind==1], okFlux$fCOS_corr[okFlux$fit_ind==1], 
     xlim=c(5,50),ylim=c(-10,30),pch=1,cex=.5)
points(okFlux$soil_temp_C_Avg[okFlux$fit_ind==2], okFlux$fCOS_corr[okFlux$fit_ind==2], 
     pch=2,cex=.5)
points(okFlux$soil_temp_C_Avg[okFlux$fit_ind==3], okFlux$fCOS_corr[okFlux$fit_ind==3], 
       pch=6,cex=.5)
points(okFlux$soil_temp_C_Avg[okFlux$fit_ind==4], okFlux$fCOS_corr[okFlux$fit_ind==4], 
       pch=5,cex=.5)
tempbins <- seq(5,50,by=0.5)
y1 <- expFunc(tempbins,coef(fit_flux_temp[[1]]))
y2 <- expFunc(tempbins,coef(fit_flux_temp[[2]]))
y3 <- expFunc(tempbins,coef(fit_flux_temp[[3]]))
y4 <- expFunc(tempbins,coef(fit_flux_temp[[4]]))
lines(tempbins, y1, lwd=2, col='darkgreen')
lines(tempbins, y2, lwd=2, col='orange')
lines(tempbins, y3, lwd=2, col='gold')
lines(tempbins, y4, lwd=2, col='red')