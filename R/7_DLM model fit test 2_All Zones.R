#ERCOT Forecasting - Basic DLM Model Fit Test on All Zones, all windows

#================================================================
# Read in load and dlm data =====================================
#================================================================

#Housekeeping
rm(list=ls())
library(matrixStats)

#Set parameters for # iterations and burn.
iter = 1100
burn = 100

#Set file path and working directory.
path = getwd()

#Set output directory.
out = paste(path,'R/Output/All Zones Fit Test/',sep='')

#Load functions.
source('R/DLM_FUNCTIONS.R')
sourceCpp(file='R/DLM_FUNCTIONS.cpp')

#Read in data.
y_all = readRDS('Data/dlm_objects/dlm_y_allzones.rda')
F_all = readRDS('Data/dlm_objects/dlm_F_allzones.rda')
G_all = readRDS('Data/dlm_objects/dlm_G_allzones.rda')

zone.names = names(F_all)

n_all = unlist(lapply(y_all,nrow))	#Sample sizes for each zone.
p_all = unlist(lapply(F_all,ncol))	#Number of predictors for each zone. 
n.zones = length(y_all)				#Number of zones.

#Save means and sds for reverse scaling later.
ybar = unlist(lapply(y_all,colMeans))
ysig = unlist(lapply(y_all,function(x) colSds(as.matrix(x)))) 

#================================================================
# Model Fit for Each Zone: ======================================
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
	y = c(scale(y))
	F = as.matrix(cbind.data.frame(int=F[,1],scale(F[,2:4]),F[,5:p]))
	
	#-----------------------------------------------------------
	#Set up hyperparameters.  (Same for all zones.)
	m0 = rep(0,p)
	C0 = diag(1,p)

	a.y = .5
	b.y = .1

	a.theta = rep(.5,p)
	b.theta = rep(.1,p)

	#-----------------------------------------------------------
	# Testing model fit for windows of data of size 10000, forecasting 100 hours into future for each window.
	zone_fit[[zone]] = dlm.fittest.movingWindow(y,F,G,K=100,m0,C0,a.y,b.y,a.theta,b.theta,win.size=10000,iter,burn)
	
	print(paste('Zone ',zone,' finished at ',Sys.time())) #Print finish time for each zone.
} 


#================================================================
#=== Plotting: ==================================================
#================================================================


#----------------------------------------------------------------
# Scaled plots
#----------------------------------------------------------------

# Loop through zones - one pdf file per zone.
for (zone in 1:n.zones){
	
	pdf(paste(out,'Scaled/',zone.names[zone],'.pdf', sep=''), height=16, width=12)

	#Set up window dimensions.
	n.win = ncol(zone_fit[[zone]]$y.pred.mean)
	par(mfrow=c( ceiling(n.win/3),3), oma=c(0,0,1.5,0), mai=c(1,.25,.25,.25))
		
	#Loop through windows.	
	for (i in 1:n.win){
		
		### Plot of forecasted values with error bars.
		t.pred 	 = seq(1,length(zone_fit[[zone]]$y.future[,i]),by=1)
		y.actual = zone_fit[[zone]]$y.future[,i]
		y.pred 	 = zone_fit[[zone]]$y.pred.mean[,i]
		lb       = y.pred - 1.96 * sqrt(zone_fit[[zone]]$y.pred.var[,i])
		ub       = y.pred + 1.96 * sqrt(zone_fit[[zone]]$y.pred.var[,i])
	
		plot(y.actual, pch=1,cex=.6, col='black', 	ylim = c(min(lb),max(ub)),
			xlab = 'time', ylab = 'load',
			main=paste('Window',i))
		
			#Shade CI region.
		polygon(c(t.pred,rev(t.pred)), c(ub, rev(lb)), col='lightgrey', border=NA)
	
		#Re-add lines and points (shaded CI covers them up).
		lines(y.pred, col='blue', lwd=2)
		lines(lb, col='lightgrey',lty=1,lwd=2)
		lines(ub, col='lightgrey',lty=1,lwd=2)
		
		points(y.actual, col='black', pch=1, cex=.6)
	} #End loop through windows plots.
	
	mtext(paste('Zone: ',zone.names[zone],sep=''), outer=T, cex=1)	
	
	dev.off()

} # End zones loop.

#----------------------------------------------------------------
# Re-Scaled plots
#----------------------------------------------------------------

# Loop through zones - one pdf file per zone.
for (zone in 1:n.zones){
	
	pdf(paste(out,'Descaled/',zone.names[zone],'.pdf', sep=''), height=16, width=12)

	#Set up window dimensions.
	n.win = ncol(zone_fit[[zone]]$y.pred.mean)
	par(mfrow=c( ceiling(n.win/3),3), oma=c(0,0,1.5,0), mai=c(1,.25,.25,.25))
		
	#Loop through windows.	
	for (i in 1:n.win){
		
		### Plot of forecasted values with error bars.
		t.pred 	 = seq(1,length(zone_fit[[zone]]$y.future[,i]),by=1)
		
		y.actual = zone_fit[[zone]]$y.future[,i] * ysig[zone] + ybar[zone]
		y.pred 	 = zone_fit[[zone]]$y.pred.mean[,i] * ysig[zone] + ybar[zone]
		
		lb       = y.pred - 1.96 * sqrt(zone_fit[[zone]]$y.pred.var[,i] * ysig[zone]^2 )
		ub       = y.pred + 1.96 * sqrt(zone_fit[[zone]]$y.pred.var[,i] * ysig[zone]^2)
	
		# ^ Since z = ay + b, where z = y scaled, b = ybar, a = 1/ysig.
		#Then var(z) = a^2 * var(y) --> var(y) = var(z) / a^2
	
		plot(y.actual, pch=1,cex=.6, col='black', 	ylim = c(min(lb),max(ub)),
			xlab = 'time', ylab = 'load',
			main=paste('Window',i))
		
			#Shade CI region.
		polygon(c(t.pred,rev(t.pred)), c(ub, rev(lb)), col='lightgrey', border=NA)
	
		#Re-add lines and points (shaded CI covers them up).
		lines(y.pred, col='blue', lwd=2)
		lines(lb, col='lightgrey',lty=1,lwd=2)
		lines(ub, col='lightgrey',lty=1,lwd=2)
		
		points(y.actual, col='black', pch=1, cex=.6)
	} #End loop through windows plots.
	
	mtext(paste('Zone: ',zone.names[zone],sep=''), outer=T, cex=1)	
	
	dev.off()

} # End zones loop.

