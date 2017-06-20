### ERCOT Forecasting - Basic DLM Model
### FUNCTIONS
### Jennifer Starling
### March 2017

library(Rcpp)
library(RcppArmadillo)

sourceCpp(file='/Users/jennstarling/UTAustin/Research/ercot/R Code/JS_DLM_FUNCTIONS.cpp')

#================================================================
# Model: ========================================================
#================================================================

#	y = F_{t} * theta_{t} + v, 		v ~ N(0,V)
#									where
#									V ~ IG(alpha_{y},beta_{y})
#									alpha_{y} = a_{y}^2 / b_{y}
#									beta_{y} = a_{y} / b_{y}
#									s.t. a_{y} = prior mean for each w_{y}
#									and b_{y} = prior variance for each w_{y}

#	theta_{t} = G * theta_{t-1}, 	w ~ N(0,W)
#									W = diag(w1...wp)

#	w1,...,wp ~ IG(alpha_i,beta_i), where
#									alpha_{i} = a_{i}^2 / b_{i}
#									beta_{i} = a_{i} / b_{i}
#									s.t. a_{i} = prior mean for each w_{i}
#									and b_{i} = prior variance for each w_{i}

#================================================================
# DLM Model Fit: ================================================
#================================================================

dlm.fit = function(y, F, F.future, G, K, m0, C0, a.y=1, b.y=2, a.theta, b.theta, iter=110, burn=10){
	#-------------------------------------------------------------
	# FUNCTION: Fits the DLM model using the gibbs sampler
	#			for ffbs using the given inputs.  Then forecasts 
	#			t+K predicted values.
	#-------------------------------------------------------------
	# INPUTS:	m0 				 = prior mean for theta.
	#			C0 				 = prior cov matrix for theta.
	#			y 				 = vector of responses
	#			F 				 = matrix of covariates, where each row Ft = vector of covariates at time t.
	#			F.future 		 = matrix of covariates for times t through K.
	#			G 				 = state matrix.
	#			a.y, b.y 		 = Hyperparams for v ~ IG(a,b). 
	#			a.theta, b.theta = Vectors of hyperparams for each w_i ~ IG(a_i,b_i).
	#			iter 			 = Number of iterations for Gibbs Sampler.
	#			burn 			 = Number of iterations to burn for Gibbs Sampler.
	#-------------------------------------------------------------	
	# OUTPUTS:	y.pred.mean		= predicted y means for t=0, t+1, t+2, ..., t+K.
	#			y.pred.var 		= predicted y vars for t=0, t+1, t+2, ..., t+K.
	#			theta.pred.mean = predicted theta vals for t=0, t+1, ..., t+K.
	#			theta.pred.var  = predicted theta cov matrices for t=0, t+1, ..., t+K.
	#			mt 				= Posterior mean of mt vector from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			Ct 				= Posterior mean of Ct matrix from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			v 				= Posterior mean of v from Gibbs Sampler.
	#			W 				= Posterior mean of W from Gibbs Sampler.
	#			theta_pm 		= Posterior mean of theta (for all times t) from Gibbs Sampler.
	#-------------------------------------------------------------
	
	#-------------------------------------------------------------
	# Data setup.
	#-------------------------------------------------------------
	
	# Number of timepoints in 'known' data set.
	T = length(y)
	
	#-------------------------------------------------------------
	# Parameter estimation of v, W, states, via Gibbs Sampler.
	#-------------------------------------------------------------
	
	# Estimate states theta_{1:T}, v and W using the Gibbs Sampler.
	# Gibbs draws v and W from full conditionals, and samples states via ffbs.
	params = gibbs(m0,C0,y,F,G,a.y,b.y,a.theta, b.theta, B=iter, burn=burn)
	
	# Save estimated parameters from Gibbs sampler.
	v = params$v_pm
	W = diag(as.numeric(params$w_pm))
	
	mt = params$mt_pm
	Ct = params$Ct_pm
	
	theta_pm = params$theta_pm
	
	#Output list for parameter estimation.
	params = list(v=v, W=W, mt=mt, Ct=Ct, theta=theta_pm)
	
	#-------------------------------------------------------------
	# In-sample fit info.
	#-------------------------------------------------------------
	
	yhat.insamp = rowSums(F * t(theta_pm))
	err.insamp  = y - yhat.insamp
	mse.insamp  = round(sum((y - yhat.insamp)^2) / length(y),8)
	
	#Output list for in-sample model fit.
	insamp = list(yhat.insamp=yhat.insamp, 
				  err.insamp=err.insamp, 
				  mse.insamp=mse.insamp)
	
	#-------------------------------------------------------------
	# Forecasting recursion.
	#-------------------------------------------------------------
	
	# Run forecasting recursion for times t+1 to t+K.
	fcast = forecast(mt,Ct,F.future,G,v,W,K)
	
	# Save forecast values.
	y.pred.mean = fcast$y_pred_mean
	y.pred.var  = fcast$y_pred_var
	
	theta.pred.mean = fcast$theta_pred_mean
	theta.pred.var  = fcast$theta_pred_var
	
	#Output list for forecasting future obs.
	forecast = list(y.pred.mean=y.pred.mean, 
					y.pred.var=y.pred.var, 
					theta.pred.mean=theta.pred.mean, 
					theta.pred.var=theta.pred.var)
		
	#-------------------------------------------------------------
	#Return output.
	#-------------------------------------------------------------
	
	return(list(params = params, 
				insamp = insamp, 
				forecast = forecast))
	
	
	# return(list(
		# y.pred.mean		= y.pred.mean,
		# y.pred.var 		= y.pred.var,
		# theta.pred.mean = theta.pred.mean,
		# theta.pred.var 	= theta.pred.var,
		# mt 				= mt,
		# Ct 				= Ct,
		# v				= v,
		# W				= W ,
		# theta_pm 		= theta_pm
	# ))	

} # End dlm.fit function.

#================================================================
# DLM Model Testing: ============================================
#================================================================

dlm.fittest = function(y, F, G, K, m0, C0, a.y=1, b.y=2, a.theta, b.theta, window, iter, burn){
	#-------------------------------------------------------------
	# FUNCTION: 	Takes entire y vector and splits it into 'known'
	#				and 'unknown' windows to predict function performance
	#				using out of sample.  Window specified by user.
	#-------------------------------------------------------------
	# INPUTS:	m0 				 = prior mean for theta.
	#			C0 				 = prior (diagonal) cov matrix for theta.
	#			y 				 = vector of responses
	#			F 				 = matrix of covariates, where each row Ft = vector of covariates at time t.
	#			F.future   		 = matrix of covariates for times t through K.
	#			G 				 = state matrix.
	#			a.y, b.y 		 = Hyperparams for v ~ IG(a,b).  a=1,b=2 = IG(.5,.5).
	#			a.theta, b.theta = Vectors of hyperparams for each w_i ~ IG(a_i,b_i).
	#			iter 			 = Number of iterations for Gibbs Sampler.
	#			burn 			 = Number of iterations to burn for Gibbs Sampler.
	#			window 			 = vector (start,stop) for what is considered the
	#							   'known' range of timepoints.  K will forecast ahead of this
	#								time range.
	#-------------------------------------------------------------	
	# OUTPUTS:	y.pred.mean 	= predicted y means for t=0, t+1, t+2, ..., t+K.
	#			y.pred.var 		= predicted y vars for t=0, t+1, t+2, ..., t+K.
	#			theta.pred.mean = predicted theta vals for t=0, t+1, ..., t+K.
	#			theta.pred.var  = predicted theta cov matrices for t=0, t+1, ..., t+K.
	#			mt 				= Posterior mean of mt vector from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			Ct 				= Posterior mean of Ct matrix from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			v 				= Posterior mean of v from Gibbs Sampler.
	#			W 				= Posterior mean of W from Gibbs Sampler.
	#			y.actual 		= actual y values for the K-step-ahead window.
	#			pred.err 		= y.pred.mean - y.actual for each predicted timepoint.
	#			MSE 			= mean squared error for range of predicted timepoints.
	#-------------------------------------------------------------

	#-------------------------------------------------------------
	# Data Setup.
	#-------------------------------------------------------------
	
	# Extract start and stop indices.
	t0 = window[1]
	tn = window[2]
	
	# Total number of time points in known data set.
	T = tn-t0 + 1	
	
	# Data setup for 'known' and 'future' observations based on window.
	y.known = y[t0:tn]
	F.known = F[t0:tn,]
	
	#First row is value at t=n+1, the time after the 'current' timepoint.
	y.future = y[(tn+1):(tn+K)]		
	F.future = F[(tn+1):(tn+K),]
	
	#-------------------------------------------------------------
	# Call DLM to predict future observations.
	#-------------------------------------------------------------
	
	fit = dlm.fit(y.known,F.known,F.future,G,K,m0,C0,a.y,b.y,a.theta,b.theta,iter,burn)

	#-------------------------------------------------------------
	# MSE Calculation.
	#-------------------------------------------------------------
	
	# Calculate MSE.pred as sum((y.pred.mean - y.future.)^2)
	y.pred = fit$forecast$y.pred.mean
	MSE.pred = mean((y.future - y.pred)^2)
	
	#-------------------------------------------------------------
	#Return output.
	#-------------------------------------------------------------
	
	params = fit$params
	insamp = fit$insamp
	forecast = fit$forecast
	
	return(list(
		params 		= params,
		insamp 		= insamp,
		forecast 	= forecast,
		y.future	= y.future,
		F.future 	= F.future,
		mse.future	= MSE.pred,
		y.known		= y.known,
		F.known		= F.known
	))

} #End dlm.fittest function.

#================================================================
# DLM Model Testing "Cross-Val" =================================
#================================================================

dlm.fittest.movingWindow = function(y, F, G, K=100, m0, C0, a.y=1, b.y=2, 
									a.theta, b.theta, win.size=10000, iter=11000, burn=1000){
	#-------------------------------------------------------------
	#FUNCTION: 	Repeats the dlm.fittest function for a moving window of 
	#			specified size.  Calculates MSE for each window, and saves
	#			predicted vs actual values for each window's fit test.
	#-------------------------------------------------------------
	#INPUTS:	m0 					= prior mean for theta.
	#			C0 					= prior (diagonal) cov matrix for theta.
	#			y 					= vector of responses
	#			F 					= matrix of covariates, where each row Ft = vector of covariates at time t.
	#			F.future 			= matrix of covariates for times t through K.
	#			G 					= state matrix.
	#			a.y, b.y 			= Hyperparams for v ~ IG(a,b).  a=1,b=2 = IG(.5,.5).
	#			a.theta, b.theta 	= Vectors of hyperparams for each w_i ~ IG(a_i,b_i).
	#			iter 				= Number of iterations for Gibbs Sampler.
	#			burn 				= Number of iterations to burn for Gibbs Sampler.
	#			win.size 			= Number of data points in each window of 'known' observations.
	#-------------------------------------------------------------	
	#OUTPUTS:	y.pred.mean = predicted y means for t=0, t+1, t+2, ..., t+K.
	#			y.actual 	= actual y values for the K-step-ahead window.
	#			pred.err 	= y.pred.mean - y.actual for each predicted timepoint.
	#			mse 		= mean squared error for range of predicted timepoints.
	#-------------------------------------------------------------
	
	#-------------------------------------------------------------
	# Windows setup.
	#-------------------------------------------------------------
	
	#Calculate number of windows. (Windows shift forward by win.size/2 each time.)
	n = length(y)

	start = 0
	end = win.size
	if(end + K > n){
		print('Error: win.size + K cannot exceed length(y).')
	}
	
	t0.vec = start	#Vector to hold start points for each window.
	tn.vec = end	#Vector to hold endpoints for each window.
	win = 1
	
	for (i in 2:n){
		start 	= t0.vec[i] = end - win.size/2
		end 	= tn.vec[i] = start + win.size
		if (end + K > n) break
		win = win + 1
	}
	
	#Trim if needed.
	t0.vec = t0.vec[1:win]
	tn.vec = tn.vec[1:win]
	
	#-------------------------------------------------------------
	#Initialize future predicted y, future actual y values, error, and mse.
	#-------------------------------------------------------------
	
	y.pred.mean = matrix(0,K,win)  #Each column is a vector of K values for a single window.
	y.pred.var  = matrix(0,K,win)  #Each column is a vector of K values for a single window.

	y.future = matrix(0,K,win)	#Each col is a vector of 'known' values for a single window.
	
	error 	= matrix(0,K,win)	#Each col is a vector of K errors for a single window.
	mse 	= rep(0,win)		#Vector of MSE values for each window.
		
	#-------------------------------------------------------------
	# Loop through each window and fit/forecast.
	#-------------------------------------------------------------
	
	#Empty list to store all dlm.fittest() objects.
	dlm.fittest.all = list()
		
	for (l in 1:win){
		
		#Run DLM calculations for window.
		output   = dlm.fittest(y, F, G, K, m0, C0, a.y, b.y, a.theta, b.theta, window=c(t0.vec[l],tn.vec[l]), iter, burn)
		
		#Save dlm fittest object.
		dlm.fittest.all[[l]] = output
		
		#Save values.
		y.pred.mean[,l] = output$forecast$y.pred.mean
		y.pred.var[,l]	= output$forecast$y.pred.var
		y.future[,l] 	= output$y.future
		error[,l] 		= output$y.future - output$forecast$y.pred.mean
		mse[l] 			= output$mse.future
	
		print(paste('Window completed:',l))
	} #End windows loop.
	
	print('Completed windows loop')
	
	#-------------------------------------------------------------
	#Return output.
	#-------------------------------------------------------------
	
	colnames(y.pred.mean) = colnames(y.pred.var) = colnames(y.future) = colnames(error) = 
		paste('window.',1:win,sep='')
	
	return(list(
		y.pred.mean		= y.pred.mean,
		y.pred.var 		= y.pred.var, 
		y.future		= y.future, 
		error			= error, 
		mse				= mse,
		dlm.fittest.all = dlm.fittest.all))	
}