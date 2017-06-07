#ERCOT Forecasting - Basic DLM Model Fit Test on Single Zone, one window.

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
sourceCpp(file='R Code/JS_DLM_FUNCTIONS.cpp')

#Set parameters for # iterations and burn.
iter = 15
burn = 10

#Read in data.
y_all = readRDS('R Data Objects/dlm_y_allzones.rda')
F_all = readRDS('R Data Objects/dlm_F_allzones.rda')
G_all = readRDS('R Data Objects/dlm_G_allzones.rda')

datetime_info = readRDS('R Data Objects/dlm_datetime_info.rda')

n_all = unlist(lapply(y_all,nrow))	#Sample sizes for each zone.
p_all = unlist(lapply(F_all,ncol))	#Number of predictors for each zone. 

n.zones = length(y_all)				#Number of zones.

ybar = unlist(lapply(y_all,colMeans))
ysig = unlist(lapply(y_all,function(x) colSds(as.matrix(x)))) 

#================================================================
# Extract Data for Zone 1: ======================================
#================================================================

### Work with Zone = COAST, Window 9:  Obs  
zone=1

### Extract zone information.
n = n_all[[zone]]
p = p_all[[zone]]

y = y_all[[zone]][,1,drop=T] 
F = F_all[[zone]]
G = G_all[[zone]]

### Select a subset of 10K 'known' plus following 100 'unknown' observations for testing.
### (COAST Window 9)

t0 = 40000
tn = 50000 - 1
K = 100

#Select subset specified above.
y = y[t0:(tn+K)]
F = F[t0:(tn+K),]
dt.info = datetime_info[t0:(tn+K),]
	
#-----------------------------------------------------------
#De-mean and scale data.
ybar = mean(y)	#Save mean for un-scaling later.
y = c(scale(y))
F = as.matrix(cbind.data.frame(int=F[,1],scale(F[,2:4]),F[,5:p]))

#-----------------------------------------------------------
# Select covariates to include.  (UPDATE COMMENT WHEN MAKING CHANGES!)
#F = F[,1:27] #Excluding daily dummies.
#p = ncol(F)
#G = diag(p)

#-----------------------------------------------------------
#Set up data subset into known/future for in-sample fit testing.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/COAST_Window 9_y.pdf')
plot(y[10000:10100],type='l', main='Coast, Window 9 (Obs 40K to 50,100)')
dev.off()

#Set 'known' and 'future' y data sets based on selected window.
t0 = 0
tn = 10000
K = 100

y.known = y[t0:tn]
F.known = F[t0:tn,]
	
y.future = y[(tn+1):(tn+K)]
F.future = F[(tn+1):(tn+K),]

length(y.known)
dim(F.known)
length(y.future)
dim(F.future)	

#================================================================
# Hyperparameter Setup: =========================================
#================================================================

#Set up hyperparameters.  (Same for all zones.)
m0 = rep(0,p)
C0 = diag(10,p)

a.y = .5
b.y = .1

a.theta = rep(.5,p)
b.theta = rep(.1,p)

### Plot density for y and theta precision priors.

alpha.y = a.y^2 / b.y 
beta.y  = a.y / b.y 

alpha.theta = a.theta^2 / b.theta
beta.theta = a.theta / b.theta

require(invgamma)

xfit = seq(.01,5,by=.001)
dens.y = dinvgamma(xfit,alpha.y,beta.y)
dens.th = dinvgamma(xfit,alpha.theta,beta.theta)

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/0_Prior_v_w_densities.pdf',height=4,width=8)
par(mfrow=c(1,2))
plot(xfit,dens.y,type='l',col='blue',lwd=2,main=paste('p(v) ~ IG(',alpha.y,',',beta.y,')',sep=''))
plot(xfit,dens.th,type='l',col='blue',lwd=2,main=paste('p(w.j) ~ IG(',alpha.theta[1],',',beta.theta[1],')',sep=''))
dev.off()

#================================================================
# In-Sample Fit Testing: ========================================
#================================================================

test.in.samp = gibbs(m0,C0,y.known,F.known,G,a.y,b.y,a.theta, b.theta, B=iter, burn)

#Calculate in-sample y values: y.t = F.t'theta.t
theta 	 = test.in.samp$theta_pm
v 		 = test.in.samp$v_pm
y.insamp = rowSums(F.known * t(theta))	#ASK ABOUT THIS FORMULA!

#In-sample MSE.
mse.insamp = round(sum((y.known - y.insamp)^2) / length(y.known),8)

#Preview results.
head(cbind(y.insamp,y.known)) 

cbind.data.frame(theta.pm.t=test.in.samp$theta_pm_t, 
	w = test.in.samp$w_pm, 
	diag.Ct = diag(test.in.samp$Ct_pm)
	)

#Plot results.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/1_Coast_Win9_Scaled_InSample.pdf',width=12,height=6)
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
### ASSESS HOURLY DUMMY VARIABLES:
### Assess potential issues with 24-hours dummy variables.
### Goal: Plot each 24-hour dummy variable's Gibbs iterations.
### Look for patterns/shark teeth where prior may be zeroing out interim obs.

#Use Ct, v, and W from gibbs sampler above.
v = test.in.samp$v_pm
W = diag(as.numeric(test.in.samp$w_pm))

Ct = test.in.samp$Ct_pm

#Run FFBS to obtain mt for all time values. (Gibbs only returns most recent.)
test.ffbs = ffbs(m0,C0,y.known,F.known,G,v,W)

#Plot each hour-of-day dummy var over time.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/2_Coast_Win9_Hour Dummies Over Time.pdf',width=18,height=12)

	par(mfrow=c(6,4),oma=c(1,1,0,0) + .1, mar=c(0,0,1,1)+1)
	for (i in 5:27){ 
		plot(test.ffbs$m[i,],type='l',xlab=paste('Hr',i-5))
		legend('topright',paste('Hr:',i-5))
	}
dev.off()

#-----------------------------------------------------------
### ASSESS DAILY DUMMY VARIABLES:
### Assess potential issues with 24-hours dummy variables.
### Goal: Plot each 24-hour dummy variable's Gibbs iterations.
### Look for patterns/shark teeth where prior may be zeroing out interim obs.

#Use Ct, v, and W from gibbs sampler above.
v = test.in.samp$v_pm
W = diag(as.numeric(test.in.samp$w_pm))

Ct = test.in.samp$Ct_pm

#Run FFBS to obtain mt for all time values. (Gibbs only returns most recent.)
test.ffbs = ffbs(m0,C0,y.known,F.known,G,v,W)

#Plot each hour-of-day dummy var over time.
pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/2_Coast_Win9_Day Dummies Over Time.pdf',width=18,height=12)

	par(mfrow=c(3,2),oma=c(1,1,0,0) + .1, mar=c(0,0,1,1)+1)
	for (i in 28:p){
		plot(test.ffbs$m[i,],type='l',xlab=paste('Day',i-27))
		legend('topright',paste('Day:',i-27))
	}
dev.off()

#================================================================
# Forecast Testing For Next 100 Obs: ============================
#================================================================

t0 = 0
tn = 10000 - 1
K=100

forecast.test = dlm.fittest(y, F, G, K, m0, C0, a.y, b.y, a.theta, b.theta, c(t0,tn), iter, burn)

#Print some results.
forecast.test$fit$mt
forecast.test$fit$Ct
forecast.test$fit$v
forecast.test$fit$W
forecast.test$fit$y.pred.mean
forecast.test$fit$y.pred.var
forecast.test$mse

#Plot results for next 100 obs.  (Predicted y means, with 95% credible interval.)

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/3_Coast_Win9_Forecast.pdf',width=12, height=6)

	par(mfrow=c(1,2))
	
	### Plot of forecasted values with error bars.
	
	t.pred 	 = seq(1,length(forecast.test$y.future),by=1)
	y.actual = forecast.test$y.future
	y.pred 	 = forecast.test$fit$y.pred.mean
	lb       = y.pred - 1.96 * sqrt(forecast.test$fit$y.pred.var)
	ub       = y.pred + 1.96 * sqrt(forecast.test$fit$y.pred.var)
	
	plot(y.actual, pch=1,cex=.6, col='black', 	ylim = c(min(lb),max(ub)),
		xlab = 'time', ylab = 'load',
		main='Coast Window 9 Fcast (Next 100 hrs)')
		
	#Shade CI region.
	polygon(c(t.pred,rev(t.pred)), c(ub, rev(lb)), col='lightgrey', border=NA)
		
	lines(y.pred, col='blue', lwd=2)
	lines(lb, col='lightgrey',lty=1,lwd=2)
	lines(ub, col='lightgrey',lty=1,lwd=2)
	
	#Re-add points (shaded CI covers them up).
	points(y.actual, col='black', pch=1, cex=.6)
	
	### Plot of residuals.
	resids = y.actual - y.pred
	qqnorm(resids,main='QQ-Norm for Forecast Errors')
	qqline(resids)
	
dev.off()

#================================================================
# PLOTS OF IN-SAMPLE ERRORS: ====================================
#================================================================

err.sq.insamp = (y.known - y.insamp)^2
day.idx = c('Mon','Tue','Wed','Thu','Fri','Sat','Sun')
day.num = 1:7
hr.num = 0:23

stat.to.plot = 'sd' # Should be 'mean' or 'sd'

cols.hrs = c('darkred','firebrick','red3','red','tomato',
	'orangered','orange3','darkorange3','darkorange',
	'darkgreen','forestgreen','green3','yellowgreen',
	'skyblue1','royalblue1','royalblue4','blue',
	'plum2','orchid','purple','purple4',
	'turquoise4 ','cyan4','darkgrey')
cols.days = c('firebrick','darkorange','limegreen','darkgreen','skyblue1','blue','purple')

#--------------------------------------------------------------
# 1: BOXPLOT: By Day for each Hour (One plot per hour)
#--------------------------------------------------------------

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/4_Coast9_InSampleMSE_Plot_1.pdf',width=30, height=40)
par(mfrow=c(6,4), oma=c(0,0,2,0))

#Loop through hours and days.  One plot per hour, over all days.
for (i in 0:23){
	
	temp.idx = dt.info$hr == i
	temp.err = err.sq.insamp[dt.info$hr == i]
	temp.day = dt.info[dt.info$hr == i,3]
	
	plot(temp.day, temp.err)
	legend('topright',paste('Hr:',i))
}

mtext("In-Sample MSE by Day for each Hour", outer = TRUE, cex = 1.5)
dev.off()

#--------------------------------------------------------------
# 2a: LINE PLOT - By Day for each Hour (One plot per hour)
#--------------------------------------------------------------

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/4_Coast9_InSampleMSE_Plot_2a.pdf',width=8, height=60)
par(mfrow=c(24,1), oma=c(0,0,2,0))

#Loop through hours and days.  One plot per hour, over all days.
for (i in 0:23){
	
	temp.idx = dt.info$hr == i
	temp.err = err.sq.insamp[temp.idx]
	temp.day = dt.info[dt.info$hr == i,3]
	temp.daynum = match(temp.day, day.idx)
	
	mse.by.day = rep(0,7)
	
	for (j in 1:7){
		
		temp.day = day.idx[j]
		temp.idx = dt.info$hr==i & dt.info$day == temp.day
		temp.err = err.sq.insamp[temp.idx]
		mse.by.day[j] = ifelse(stat.to.plot=='mean',mean(temp.err,na.rm=T),sd(temp.err,na.rm=T))
	}
	
	plot(1:7, mse.by.day, type='l', col='blue', lwd=2, ylab='In Sample MSE', xlab='Day', xaxt='n')
	axis(1, at=1:7, labels=day.idx)
	legend('topright',paste('Hr:',i))
}

mtext(paste("In-Sample MSE by Day for each Hour:",stat.to.plot), outer = TRUE, cex = 1.5)
dev.off()

#--------------------------------------------------------------
# 2b: LINE PLOTS - By Day for each Hour, Overlaid On One Plot
#--------------------------------------------------------------

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/4_Coast9_InSampleMSE_Plot_2b.pdf',width=14, height=9)
par(mfrow=c(1,1), oma=c(0,0,2,0), xpd=T, mar = c(5.1, 4.1, 4.1, 8.1))

#Loop through hours and days.  One plot per hour, over all days.
for (i in 0:23){
	
	temp.idx = dt.info$hr == i
	temp.err = err.sq.insamp[temp.idx]
	temp.day = dt.info[dt.info$hr == i,3]
	temp.daynum = match(temp.day, day.idx)
	
	mse.by.day = rep(0,7)
	
	for (j in 1:7){
		
		temp.day = day.idx[j]
		temp.idx = dt.info$hr==i & dt.info$day == temp.day
		temp.err = err.sq.insamp[temp.idx]
		mse.by.day[j] = ifelse(stat.to.plot=='mean',mean(temp.err,na.rm=T),sd(temp.err,na.rm=T))
	}
	
	if(i==0){
		plot(1:7, mse.by.day, type='l', col=cols.hrs[1], lwd=3, ylab='In Sample MSE', xlab='Day', xaxt='n', 
			ylim=c(0,.008))
		axis(1, at=1:7, labels=day.idx)
		legend('topright',as.character(hr.num), lty=rep(1,24), lwd=rep(4,24), col=cols.hrs, inset=c(-0.1,0), title='Hour')
	} else{
		lines(1:7, mse.by.day, col=cols.hrs[i+1], lwd=3)	
	}
}

mtext(paste("In-Sample MSE by Day for each Hour:",stat.to.plot), outer = TRUE, cex = 1.5)
dev.off()

#--------------------------------------------------------------
# 3a: LINE PLOTS - By Hour for each Day
#--------------------------------------------------------------

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/Coast9_InSampleMSE_Plot_3a.pdf',width=8, height=20)
par(mfrow=c(7,1), oma=c(0,0,2,0))

#Loop through hours and days.  One plot per hour, over all days.
for (j in 1:7){
	
	temp.day = day.idx[j]
	mse.by.hr = rep(0,24)
	
	for (i in 1:24){
		temp.hr = i-1
		temp.idx = dt.info$day==temp.day & dt.info$hr==temp.hr
		temp.err = err.sq.insamp[temp.idx]
		mse.by.hr[i] = ifelse(stat.to.plot=='mean',mean(temp.err,na.rm=T),sd(temp.err,na.rm=T))
	}
	
	plot(1:24, mse.by.hr, type='l', col='blue', lwd=2, ylab='In Sample MSE', xlab='Hour', xaxt='n')
	axis(1, at=1:24, labels=0:23)
	legend('topright',paste('Day:',temp.day))
}

mtext(paste("In-Sample MSE by Hour for Each Day:",stat.to.plot), outer = TRUE, cex = 1.5)
dev.off()

#--------------------------------------------------------------
# 3b: LINE PLOTS - By Hour for each Day, Overlaid On One Plot
#--------------------------------------------------------------

pdf('/Users/jennstarling/UTAustin/Research/ercot/Figures/Small Scale Fit Test/Coast9_InSampleMSE_Plot_3b.pdf',width=14, height=9)
par(mfrow=c(1,1), oma=c(0,0,2,0), xpd=T, mar = c(5.1, 4.1, 4.1, 8.1))

#Loop through hours and days.  One plot per hour, over all days.
for (j in 1:7){
	
	temp.day = day.idx[j]
	mse.by.hr = rep(0,24)
	
	for (i in 1:24){
		temp.hr = i-1
		temp.idx = dt.info$day==temp.day & dt.info$hr==temp.hr
		temp.err = err.sq.insamp[temp.idx]
		mse.by.hr[i] = ifelse(stat.to.plot=='mean',mean(temp.err,na.rm=T),sd(temp.err,na.rm=T))
	}
	
	if(j==1){
		plot(1:24, mse.by.hr, type='l', col=cols.days[j], lwd=2, ylab='In Sample MSE', xlab='Hour', xaxt='n')
		axis(1, at=1:24, labels=0:23)
		legend('topright',day.idx, lty=rep(1,7), lwd=rep(4,7), col=cols.days, inset=c(-0.1,0), title='Day')
	} else{
		lines(1:24, mse.by.hr, col=cols.days[j], lwd=2)
	}
}

mtext(paste("In-Sample MSE by Hour for Each Day:",stat.to.plot), outer = TRUE, cex = 1.5)
dev.off()
