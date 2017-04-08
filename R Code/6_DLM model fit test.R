#ERCOT Forecasting - Basic DLM Model

#================================================================
# Read in load and dlm data =====================================
#================================================================

#Housekeeping
rm(list=ls())
library(matrixStats)

#Set working directory.
setwd('/Users/jennstarling/UTAustin/Research/ercot')

#Load functions.
source('R Code/JS_DLM_FUNCTIONS.R')
sourceCpp(file='R Code/JS_DLM_FUNCTIONS_2.cpp')

#Read in data.
y_all = readRDS('R Data Objects/dlm_y_allzones.rda')
F_all = readRDS('R Data Objects/dlm_F_allzones.rda')
G_all = readRDS('R Data Objects/dlm_G_allzones.rda')

n_all = unlist(lapply(y_all,nrow))	#Sample sizes for each zone.
p_all = unlist(lapply(F_all,ncol))	#Number of predictors for each zone. 

n.zones = length(y_all)				#Number of zones.

ybar = unlist(lapply(y_all,colMeans))
ysig = unlist(lapply(y_all,function(x) colSds(as.matrix(x)))) 

#================================================================
# Check Model Fit: ==============================================
#================================================================

#Empty list to hold fit info for each zone.
zone_fit = list()

Sys.time()	#Begin system time.

#Loop through zones.
for (zone in 1:n.zones){
	
	### Extract zone information.

	#Extract DLM known data for selected zone.
	n = n_all[[zone]]
	p = p_all[[zone]]

	y = y_all[[zone]][,1,drop=T] 
	F = F_all[[zone]]
	G = G_all[[zone]]
	
	#-----------------------------------------------------------
	#De-mean and scale data.
	ybar = mean(y)	#Save mean for un-scaling later.
	y = c(scale(y))
	F = as.matrix(cbind.data.frame(int=F[,1],scale(F[,2:4]),F[,5:p]))
	
	#-----------------------------------------------------------
	#Set up hyperparameters.  (Same for all zones.)
	m0 = rep(0,p)
	C0 = diag(10,p)

	# v ~ IG(.5,.5)
	a.y = 1
	b.y = 2	

	# w_i ~ IG(.5,.5)
	a.theta = rep(1,p)
	b.theta = rep(2,p)

	#-----------------------------------------------------------
	# Testing model fit for windows of data of size 10000, forecasting 100 hours into future for each window.
	zone_fit[[zone]] = dlm.fittest.cv(y,F,G,K=100,m0,C0,a.y=1,b.y=2,a.theta,b.theta,win.size=10000,iter=1100,burn=1000)
	
	print(paste('Zone ',zone,' finished at ',Sys.time())) #Print finish time for each zone.
} 

#COAST PLOTTING
zone_fit[[1]]$mse
mean(zone_fit[[1]]$mse)

cbind(zone_fit[[1]]$y.pred,zone_fit[[1]]$y.known)

#Plot un-rescaled data.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Coast_Scaled.pdf')
par(mfrow=c(3,3)) #One plot per window.
for (i in 1:9){
	y.p = zone_fit[[1]]$y.pred[,i]
	y.k = zone_fit[[1]]$y.known[,i]
	err = zone_fit[[1]]$error[,i]
	
	plot(y.p,col='blue',type='l',ylim=c(-3,3),main=paste('Coast, window',i))
	points(y.k,col='black',type='l')
}
dev.off()

#Rescaled plots.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Coast_UnScaled.pdf')
par(mfrow=c(3,3)) #One plot per window.
for (i in 1:9){
	yp.rescaled = zone_fit[[1]]$y.pred[,i]  * ysig[[1]] + ybar[[1]]
	yk.rescaled = zone_fit[[1]]$y.known[,i] * ysig[[1]] + ybar[[1]]
	err = yp.rescaled - yk.rescaled
	mse = mean(err^2)
	
	print(mse)
	plot(yp.rescaled,col='blue',type='l',ylim=c(0,30000),main=paste('Coast, window',i))
	points(yk.rescaled,col='black',type='l')
}
dev.off()

#================================================================
# Fit Troubleshooting - Small Scale: ============================
#================================================================

#Set parameters for # iterations and burn.
iter = 1100
burn = 100

### Work with Zone = COAST, Window 9:  Obs  
zone=1

### Select a subset of 10K 'known' plus following 100 'unknown' observations for testing.
### (COAST Window 9)

t0 = 40000
tn = 50000 - 1
K = 100

#-----------------------------------------------------------
### Extract zone information.

	#Extract DLM known data for selected zone.
	n = n_all[[zone]]
	p = p_all[[zone]]

	y = y_all[[zone]][,1,drop=T] 
	F = F_all[[zone]]
	G = G_all[[zone]]
	
	#Select subset specified above.
	y = y[t0:(tn+K)]
	F = F[t0:(tn+K),]
	
	#-----------------------------------------------------------
	#De-mean and scale data.
	ybar = mean(y)	#Save mean for un-scaling later.
	y = c(scale(y))
	F = as.matrix(cbind.data.frame(int=F[,1],scale(F[,2:4]),F[,5:p]))
	
	#-----------------------------------------------------------
	#Plot y values to be forecasted.
	plot(y[10000:10100],type='l')
		
	#-----------------------------------------------------------
	#Set up hyperparameters.  (Same for all zones.)
	m0 = rep(0,p)
	C0 = diag(1,p)

	# v ~ IG(.5,.5)
	a.y = 1
	b.y = 2	

	# w_i ~ IG(.5,.5)
	a.theta = rep(1,p)
	b.theta = rep(2,p)

#-----------------------------------------------------------
### ASSESS HOURLY DUMMY VARIABLES:
### Assess potential issues with 24-hours dummy variables.
### Goal: Plot each 24-hour dummy variable's Gibbs iterations.
### Look for patterns/shark teeth where prior may be zeroing out interim obs.

#First run Gibbs Sampler to obtain v and W.
test.gibbs = gibbs(m0,C0,y.known,F.known,G,a.y,b.y,a.theta, b.theta, B=iter, burn)
v = test.gibbs$v_pm
W = diag(as.numeric(test.gibbs$w_pm))

#Run FFBS to obtain mt for all time values. (Gibbs only returns most recent.)
test.ffbs = ffbs(m0,C0,y.known,F.known,G,v,W,bs=1)

#Plot each hour-of-day dummy var over time.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Coast_Win9_Hour Dummies Over Time.pdf',width=18,height=12)
par(mfrow=c(6,4),oma=c(1,1,0,0) + .1, mar=c(0,0,1,1)+1)
for (i in 1:24){
	plot(test.ffbs$m_all_t[i,],type='l',xlab=paste('Hr',i-1))
	legend('topright',paste('Hr:',i-1))
}
dev.off()

#-----------------------------------------------------------
### IN-SAMPLE FIT:
### Try an in-sample fit plot using the same v and W as above.

#Calculate in-sample y values: y.t = F.t'theta.t
theta = test.gibbs$theta_pm
y.insamp = rowSums(F.known * t(theta)) 	

#In-sample MSE.
mse.insamp = round(sum((y.known - y.insamp)^2) / length(y.known),8)

#Preview results.
head(cbind(y.insamp,y.known)) 

#Plot results.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Coast_Win9_Scaled_InSample.pdf',width=12,height=6)
	par(mfrow=c(1,2))

	#Plot fit.
	plot(y.known[seq(1,10000,by=100)],type='l',
		main = paste('Coast Win 9, In-Sample, MSE = ',mse.insamp))
	lines(y.insamp[seq(1,10000,by=100)],col='blue')
	legend("top",c(y.known,y.insamp),c('y','y.pred'),col=c('black','blue'),lty=c(1,1))

	#Plot residuals.
	err = y.known - y.insamp
	qqnorm(err,main='QQ-Norm for In-Sample Errors')
	qqline(err)
dev.off()

#-----------------------------------------------------------
### PLOT PREDICTED VALUES:
### Try plotting predicted means and predicted draws for y.

# Data setup.
t0 = 0
tn = 10000
K = 100

y.known = y[t0:tn]
F.known = F[t0:tn,]
	
y.future = y[(tn+1):(tn+K)]
F.future = F[(tn+1):(tn+K),]

# Run DLM fit and save predicted values.
coast_win9 = dlm.fit(y.known,F.known,F.future,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,iter,burn)
y.pr = coast_win9$y.pred.draw

#v: 0.000156267
#range(w): 0.000000e+00, 1.367961e-08

# Plot results.

#Plot predicted means.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Coast_Win9_Scaled_Predicted.pdf',width=12,height=4)
par(mfrow=c(1,3))
	
	#Plot actual future y values vs mean predicted y values, with CI.
	yp.m = coast_win9$y.pred.mean
	yp.v = coast_win9$y.pred.var
	lb = yp.m - 1.96*sqrt(yp.v)
	ub = yp.m + 1.96*sqrt(yp.v)
	plot(yp.m,type='l',col='blue',ylim=range(min(lb),max(ub)),xlab='time',ylab='load',
		main='Predicted Future Y Means vs Observed Future Y')
	lines(y.future)
	lines(lb,lty=2,col='red')
	lines(ub,lty=2,col='red')
	legend("topright",c(y.future,yp.m,lb),c('y','y.pred.mean','error bars'),col=c('black','blue','red'),lty=c(1,1,2),bg = 'white')
	
	#Plot actual future y values vs predicted y draws.
	yp.d = coast_win9$y.pred.draw
	plot(yp.d,type='l',col='blue',main = 'Predicted Future Y Draws vs Observed Future Y')
	lines(y.future)
	legend("top",c(y.future,yp.d),c('y','y.pred.draw'),col=c('black','blue'),lty=c(1,1),bg = 'white')
	
	# Plot residuals; check for normality.
	errs = y.future - yp.d
	qqnorm(errs,main='QQ-Norm: Forecast Err (Y Draw vs Observed Future Y)')
	qqline(errs)
dev.off()
#-----------------------------------------------------------

