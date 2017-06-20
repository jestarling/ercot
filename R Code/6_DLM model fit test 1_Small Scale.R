#ERCOT Forecasting - Basic DLM Model Fit Test on Single Zone, one window.

#=======================================================================
#===   Read in load and dlm data   =====================================
#=======================================================================

#Housekeeping
rm(list=ls())
library(matrixStats)

#Save script start time.
time.start = Sys.time()

#Set file path and working directory.
path = '/Users/jennstarling/UTAustin/Research/ercot'
setwd(path)

#Set output directory.
out = paste(path,'/Figures/Small Scale Fit Test/',sep='')

#Load functions.
source('R Code/JS_DLM_FUNCTIONS.R')
sourceCpp(file='R Code/JS_DLM_FUNCTIONS.cpp')

#Set parameters for # iterations and burn.
iter = 1100
burn = 100

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

#=======================================================================
#=== Extract Data for Zone 1: ==========================================
#=======================================================================

#-----------------------------------------------------------
#   Select data from Zone = COAST, Window 9:  Obs  
#-----------------------------------------------------------

zone=1

### Extract zone information.
n = n_all[[zone]]
p = p_all[[zone]]

y = y_all[[zone]][,1,drop=T] 
F = F_all[[zone]]
G = G_all[[zone]]

### Select a window subset of 10K 'known' plus following 100 'unknown' 
### observations for testing. (COAST Window 9)

t0 = 40000
tn = 50000 - 1
K = 100

#Select subset specified above as the 'full' data set of known and future.
y = y[t0:(tn+K)]
F = F[t0:(tn+K),]
dt.info = datetime_info[t0:(tn+K),]
	
#-----------------------------------------------------------
#   De-mean and scale data.
#-----------------------------------------------------------

ybar = mean(y)	#Save mean for un-scaling later.
ysig = sd(y)	#Save var for un-scaling later.

y = c(scale(y))
F = as.matrix(cbind.data.frame(int=F[,1],scale(F[,2:4]),F[,5:p]))

#-----------------------------------------------------------
#   Set up data subset into known/future for in-sample fit testing.
#-----------------------------------------------------------

#Set 'known' and 'future' y data sets based on selected window.
t0 = 0
tn = 10000
win = c(t0,tn)
K = 100

#Save y, F and date/time info for each obs for 'known' data.
y.known = y[t0:tn]
F.known = F[t0:tn,]
dt.info.known = dt.info[t0:tn,]
	
#Save y, F and date/time info for each obs for 'future' data.	
y.future = y[(tn+1):(tn+K)]
F.future = F[(tn+1):(tn+K),]
dt.info.future = dt.info[(tn+1):(tn+K),]

#Confirm dimensions are as expected.
length(y.known)
dim(F.known)
length(y.future)
dim(F.future)	

#Plot 'known' and 'future' data.
pdf(paste(out,"01_COAST_Window 9_Raw Data (y).pdf"), width=12, height=6)
	par(mfrow=c(1,2))
	
	plot(y.known, type='l', xlab='t', main='Coast, Window 9, Known (40K to 50K)')
	plot(y.future, type='l', xlab='t', main='Coast, Window 9, Future (50K to 50,100)')
	
dev.off()

#=======================================================================
#===   Hyperparameter Setup:   =========================================
#=======================================================================

#Set up hyperparameters.  (Same for all zones.)
m0 = rep(0,p)
C0 = diag(1,p)

a.y = .5
b.y = .1

a.theta = rep(.5,p)
b.theta = rep(.1,p)

#-----------------------------------------------------------
#   Plot density for y and theta precision priors.
#-----------------------------------------------------------

alpha.y = a.y^2 / b.y 
beta.y  = a.y / b.y 

alpha.theta = a.theta^2 / b.theta
beta.theta = a.theta / b.theta

require(invgamma)

xfit = seq(.01,5,by=.001)
dens.y = dinvgamma(xfit,alpha.y,beta.y)
dens.th = dinvgamma(xfit,alpha.theta,beta.theta)

pdf(paste(out,"02_Prior_v_w_densities.pdf"), height=4, width=8)
	par(mfrow=c(1,2))
	plot(xfit,dens.y,type='l',col='blue',lwd=2,main=paste('p(v) ~ IG(',alpha.y,',',beta.y,')',sep=''))
	plot(xfit,dens.th,type='l',col='blue',lwd=2,main=paste('p(w.j) ~ IG(',alpha.theta[1],',',beta.theta[1],')',sep=''))
dev.off()

#=======================================================================
#== Run DLMFIT for in-sample and forecast:   ===========================
#=======================================================================

output = dlm.fittest(y, F, G, K, m0, C0, a.y, b.y, a.theta, b.theta, window=c(t0,tn), iter, burn)

#=======================================================================
#== In-Sample Fit Testing: =============================================
#=======================================================================

#-----------------------------------------------------------
#   Extract, preview and save parameters from Gibbs Sampler.
#-----------------------------------------------------------

v = output$params$v
W = output$params$W
w = diag(W)
mt = output$params$mt
Ct = output$params$Ct
theta = output$params$theta
theta.t = theta[,ncol(theta)]

param.post.means = cbind.data.frame(
	theta.pm.t = theta.t,
	w = w,
	diag.Ct = diag(Ct))
	
write.csv(param.post.means, paste(out,'params_pm.csv'))	

#----------------------------------------------------
#   Fit model to data, in-sample.
#----------------------------------------------------

#Extract results.
y.insamp 	= output$insamp$yhat.insamp
mse.insamp 	= output$insamp$mse.insamp
err.insamp  = output$insamp$err.insamp

#Preview results.
mse.insamp
head(cbind(y.insamp,y.known)) 

#----------------------------------------------------
#   Plot results of in-sample fit.
#----------------------------------------------------

pdf(paste(out,"03_COAST_Win9_Scaled_InSampleFit.pdf"), width=12,height=6)

	par(mfrow=c(1,2))

	#Plot fit.
	plot(y.known[seq(1,10000,by=100)],type='l',
		main = paste('Coast Win 9, In-Sample, MSE = ',mse.insamp))
	lines(y.insamp[seq(1,10000,by=100)],col='blue')
	legend("top",c(y.known,y.insamp),c('y','y.pred'),col=c('black','blue'),lty=c(1,1))

	#Plot residuals.
	qqnorm(err.insamp,main='QQ-Norm for In-Sample Errors')
	qqline(err.insamp)
dev.off()

#----------------------------------------------------
#   Trace Plots: Temp, Temp^2, Bushr 
#	Plot each variable's trace over time.       
#----------------------------------------------------

#Run FFBS to obtain mt for all time values. (Gibbs only returns most recent.)
test.ffbs = ffbs(m0,C0,y.known,F.known,G,v,W)

#Plot each hour-of-day dummy var over time.
pdf(paste(out,"05_COAST_Win9_Bushr and Temp Trace.pdf"), width=18, height=12)
	
	bustemp.cols = c(which(substr(colnames(F),3,6)=='temp'), which( colnames(F)=="F.holiday"))
	bustemp.colnames = colnames(F)[bustemp.cols]
	par(mfrow=c(2,2),oma=c(1,1,0,0) + .1, mar=c(0,0,1,1)+1)
	
	for (i in 1:length(bustemp.cols)){ 
		plot(test.ffbs$m[bustemp.cols[i],],type='l',xlab=bustemp.colnames[i])
		legend('topright',bustemp.colnames[i])
	}
dev.off()

#--------------------------------------------------------
#   Trace Plots: Hourly Dummies
#	Plot each variable's trace over time.   
#--------------------------------------------------------

#Run FFBS to obtain mt for all time values. (Gibbs only returns most recent.)
test.ffbs = ffbs(m0,C0,y.known,F.known,G,v,W)

#Plot each hour-of-day dummy var over time.
pdf(paste(out,"06_COAST_Win9_Hourly Dummy Trace.pdf"), width=18, height=12)
	
	hr.cols = which(substr(colnames(F),1,2)=="hr")
	par(mfrow=c(6,4),oma=c(1,1,0,0) + .1, mar=c(0,0,1,1)+1)
	
	for (i in hr.cols){ 
		plot(test.ffbs$m[i,],type='l',xlab=paste('Hr',i-hr.cols[1]))
		legend('topright',paste('Hr:',i-hr.cols[1]))
	}
dev.off()

#--------------------------------------------------------
#   Plot autocorrelation of in-sample errors. 
#--------------------------------------------------------

#Save in-sample errors.
err.insamp = output$insamp$err.insamp

#Confirm errors are mean zero.
mean(err.insamp)

#Plot in-sample errors' autocorrelation structure.
pdf(paste(out,"Insamp_Errors_Autocorrelation.pdf"), width=10, height=8)
	acf(err.insamp)
dev.off()

#=======================================================================
#===   Forecast Testing For Next 100 Obs:   ============================
#=======================================================================

# Uses forecast output from dlm.fittest run.

#--------------------------------------------------------
#   Extract forecast results.
#--------------------------------------------------------

y.pred.mean = output$forecast$y.pred.mean
y.pred.var  = output$forecast$y.pred.var

theta.pred.mean = output$forecast$theta.pred.mean
theta.pred.var = output$forecast$theta.pred.var

y.future = output$y.future
mse.pred = output$mse.future

#--------------------------------------------------------
#   Plot results for next 100 obs.  
#   (Predicted y means, with 95% credible interval.)
#--------------------------------------------------------

pdf(paste(out,"04a_COAST_Win9_Forecast.pdf"), width=12, height=6)

	par(mfrow=c(1,2))
	
	### Plot of forecasted values with error bars.
	
	t.pred 	 = seq(1,length(y.future),by=1)
	
	lb = y.pred.mean - 1.96 * sqrt(y.pred.var)
	ub = y.pred.mean + 1.96 * sqrt(y.pred.var)
	
	plot(y.future, pch=1,cex=.6, col='black', ylim = c(min(lb),max(ub)),
		xlab = 'time', ylab = 'load',
		main='Coast Window 9 Fcast (Next 100 hrs)')
		
	#Shade CI region.
	polygon(c(t.pred,rev(t.pred)), c(ub, rev(lb)), col='lightgrey', border=NA)
		
	lines(y.pred.mean, col='blue', lwd=2)
	lines(lb, col='lightgrey',lty=1,lwd=2)
	lines(ub, col='lightgrey',lty=1,lwd=2)
	
	#Re-add points (shaded CI covers them up).
	points(y.future, col='black', pch=1, cex=.6)
	
	### Plot of residuals.
	resids = y.future - y.pred.mean
	qqnorm(resids,main='QQ-Norm for Forecast Errors')
	qqline(resids)
	
dev.off()

#--------------------------------------------------------
#   Plot results in original scale (rescaled) for next 100 obs.  
#   (Predicted y means, with 95% credible interval.)
#--------------------------------------------------------

y.pred.mean.rescaled = y.pred.mean * ysig + ybar
y.pred.var.rescaled  = y.pred.var * ysig^2 
	#Since z = ay + b, where z = y scaled, b = ybar, a = 1/ysig.
	#Then var(z) = a^2 * var(y) --> var(y) = var(z) / a^2

theta.pred.mean = output$forecast$theta.pred.mean
theta.pred.var = output$forecast$theta.pred.var

y.future.rescaled = y.future*ysig + ybar
mse.pred.rescaled = mean((y.pred.mean.rescaled - y.future.rescaled)^2)

pdf(paste(out,"04b_COAST_Win9_Forecast_Rescaled.pdf"), width=12, height=6)

	par(mfrow=c(1,2))
	
	### Plot of forecasted values with error bars.
	
	t.pred 	 = seq(1,length(y.future),by=1)
	
	lb = y.pred.mean.rescaled - 1.96 * sqrt(y.pred.var.rescaled)
	ub = y.pred.mean.rescaled + 1.96 * sqrt(y.pred.var.rescaled)
	
	plot(y.future.rescaled, pch=1,cex=.6, col='black', ylim = c(min(lb),max(ub)),
		xlab = 'time', ylab = 'load',
		main='Coast Window 9 Fcast (Next 100 hrs)')
		
	#Shade CI region.
	polygon(c(t.pred,rev(t.pred)), c(ub, rev(lb)), col='lightgrey', border=NA)
		
	#Re-add points and lines (shaded CI covers them up).
	lines(lb, col='lightgrey',lty=1,lwd=2)
	lines(ub, col='lightgrey',lty=1,lwd=2)
	points(y.future.rescaled, col='black', pch=1, cex=.6)
	lines(y.pred.mean.rescaled, col='blue', lwd=2)
	
	### Plot of residuals.
	resids = y.future.rescaled - y.pred.mean.rescaled
	qqnorm(resids,main='QQ-Norm for Forecast Errors')
	qqline(resids)
	
dev.off()


#=======================================================================
#===   PLOTS OF IN-SAMPLE MSE:   =======================================
#=======================================================================

err.sq.insamp = (y.known - y.insamp)^2
day.idx = c('Mon','Tue','Wed','Thu','Fri','Sat','Sun')
day.num = 1:7
hr.num = 0:23

cols.hrs = c('darkred','firebrick','red3','red','tomato',
	'orangered','orange3','darkorange3','darkorange',
	'darkgreen','forestgreen','green3','yellowgreen',
	'skyblue1','royalblue1','royalblue4','blue',
	'plum2','orchid','purple','purple4',
	'turquoise4 ','cyan4','darkgrey')
cols.days = c('firebrick','darkorange','limegreen','darkgreen','skyblue1','blue','purple')

#--------------------------------------------------------
#   1: BOXPLOT: By Day for each Hour (One plot per hour)
#--------------------------------------------------------

pdf(paste(out,"InsampMSE_01_COAST_Win9_Boxplot By Day for ea Hr.pdf"), width=15, height=25)

	par(mfrow=c(6,4), oma=c(0,0,2,0))

	#Loop through hours and days.  One plot per hour, over all days.
	for (i in 0:23){
	
		temp.idx = dt.info$hr == i
		temp.err = err.sq.insamp[dt.info$hr == i]
		temp.day = dt.info[dt.info$hr == i,3]
	
		plot(temp.day, temp.err, ylim = c(0,quantile(err.sq.insamp,.975)) )
		legend('topright',paste('Hr:',i))
	}

	mtext("In-Sample MSE by Day for each Hour", outer = TRUE, cex = 1.5)
dev.off()

#--------------------------------------------------------
#   2: BOXPLOT: By Hour for each Day (One plot per day)
#--------------------------------------------------------

pdf(paste(out,"InsampMSE_02_COAST_Win9_Boxplot By Hr for ea Day.pdf"), width=15, height=25)

	par(mfrow=c(7,1), oma=c(0,0,2,0))
	days = c('Mon','Tue','Wed','Thu','Fri','Sat','Sun')

	#Loop through hours and days.  One plot per hour, over all days.
	for (i in 1:7){
	
		temp.idx = dt.info$day == days[i]
		temp.err = err.sq.insamp[temp.idx]
		temp.hr = dt.info[temp.idx,2]
	
		#plot(temp.hr, temp.err)
		boxplot(temp.err ~ temp.hr, ylim = c(0,quantile(err.sq.insamp,.975)) )
		legend('topright',paste('Day: ',days[i]))
	}

	mtext("In-Sample MSE by Hour for each Day", outer = TRUE, cex = 1.5)
dev.off()

#=======================================================================
#===   OUTPUT PERFORMANCE METRICS:   ===================================
#=======================================================================

#Save script end time.
time.end = Sys.time()
time.run = time.end - time.start

#Print info on run-time and iterations.
time.start
time.end
time.run
iter

system("say DLM script number 1 has completed!")
