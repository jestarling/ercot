### ERCOT Forecasting - Basic DLM Model
### FUNCTIONS
### Jennifer Starling
### March 2017

library(Rcpp)
library(RcppArmadillo)

sourceCpp(file='/Users/jennstarling/UTAustin/Research/ercot/R Code/JS_DLM_FUNCTIONS_2.cpp')

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

dlm.fit = function(y,F,F.future,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,iter=11000,burn=1000){
	#-------------------------------------------------------------
	#FUNCTION: 	Fits the DLM model using the gibbs sampler
	#			for ffbs using the given inputs.  Then forecasts 
	#			t+K predicted values.
	#-------------------------------------------------------------
	#INPUTS:	m0 = prior mean for theta.
	#			C0 = prior (diagonal) cov matrix for theta.
	#			y = vector of responses
	#			F = matrix of covariates, where each row
	#				Ft = vector of covariates at time t.
	#			F.future = matrix of covariates for times t through K.
	#			G = state matrix.
	#			a.y, b.y = Hyperparams for v ~ IG(a,b).  a=1,b=2 = IG(.5,.5).
	#			a.theta, b.theta = Vectors of hyperparams for each w_i ~ IG(a_i,b_i).
	#			iter = Number of iterations for Gibbs Sampler.
	#			burn = Number of iterations to burn for Gibbs Sampler.
	#-------------------------------------------------------------	
	#OUTPUTS:	y.pred.mean = predicted y means for t=0, t+1, t+2, ..., t+K.
	#			y.pred.var = predicted y vars for t=0, t+1, t+2, ..., t+K.
	#			theta.pred.mean = predicted theta vals for t=0, t+1, ..., t+K.
	#			theta.pred.var = predicted theta cov matrices for t=0, t+1, ..., t+K.
	#			mt = Posterior mean of mt vector from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			Ct = Posterior mean of Ct matrix from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			v = Posterior mean of v from Gibbs Sampler.
	#			W = Posterior mean of W from Gibbs Sampler.
	#-------------------------------------------------------------
	
	#Estimate parameters using ffbs via Gibbs Sampler.
	params = gibbs(m0,C0,y,F,G,a.y,b.y,a.theta, b.theta, B=iter, burn=burn)
	
	#Forecast using estimated params.  Forecast includes current time t+1 to t+K.
	
	mt = params$mt_pm
	Ct = params$Ct_pm
	v = params$v_pm
	W = diag(as.numeric(params$w_pm))
	
	fcast = forecast(mt,Ct,F.future,G,v,W,K)
	
	y.pred.draw = fcast$y_pred_draw
	y.pred.mean = fcast$y_pred_mean
	y.pred.var = fcast$y_pred_var
	theta.pred.mean = fcast$theta_pred_mean
	theta.pred.var = fcast$theta_pred_var
	
	return(list(
		y.pred.draw=y.pred.draw,
		y.pred.mean=y.pred.mean,
		y.pred.var = y.pred.var,
		theta.pred.mean = theta.pred.mean,
		theta.pred.var = theta.pred.var,
		mt = mt,
		Ct = Ct,
		v=v,
		W=W))	
}

#================================================================
# DLM Model Testing: ============================================
#================================================================

dlm.fittest = function(y,F,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,window,iter,burn){
	#-------------------------------------------------------------
	#FUNCTION: 	Takes entire y vector and splits it into 'known'
	#			and 'unknown' windows to predict function performance
	#			using out of sample.  Window specified by user.
	#-------------------------------------------------------------
	#INPUTS:	m0 = prior mean for theta.
	#			C0 = prior (diagonal) cov matrix for theta.
	#			y = vector of responses
	#			F = matrix of covariates, where each row
	#				Ft = vector of covariates at time t.
	#			F.future = matrix of covariates for times t through K.
	#			G = state matrix.
	#			a.y, b.y = Hyperparams for v ~ IG(a,b).  a=1,b=2 = IG(.5,.5).
	#			a.theta, b.theta = Vectors of hyperparams for each w_i ~ IG(a_i,b_i).
	#			iter = Number of iterations for Gibbs Sampler.
	#			burn = Number of iterations to burn for Gibbs Sampler.
	#			window = vector (start,stop) for what is considered the
	#				'known' range of timepoints.  K will forecast ahead of this
	#				time range.
	#-------------------------------------------------------------	
	#OUTPUTS:	y.pred.mean = predicted y means for t=0, t+1, t+2, ..., t+K.
	#			y.pred.var = predicted y vars for t=0, t+1, t+2, ..., t+K.
	#			theta.pred.mean = predicted theta vals for t=0, t+1, ..., t+K.
	#			theta.pred.var = predicted theta cov matrices for t=0, t+1, ..., t+K.
	#			mt = Posterior mean of mt vector from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			Ct = Posterior mean of Ct matrix from Gibbs Sampler. (Mean of ffbs iterations for time t.)
	#			v = Posterior mean of v from Gibbs Sampler.
	#			W = Posterior mean of W from Gibbs Sampler.
	#			y.actual = actual y values for the K-step-ahead window.
	#			pred.err = y.pred.mean - y.actual for each predicted timepoint.
	#			MSE = mean squared error for range of predicted timepoints.
	#-------------------------------------------------------------

	#Extract start and stop indices.
	t0 = window[1]
	tn = window[2]
	T = tn-t0 + 1	#Total number of time points in known data set.
	
	y.known = y[t0:tn]
	F.known = F[t0:tn,]
	
	y.future = y[(tn+1):(tn+K)]		#First row is value at t=n, the 'current' timepoint.
	F.future = F[(tn+1):(tn+K),]
	
	fit = dlm.fit(y.known,F.known,F.future,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,iter,burn)
	
	#Don't want t=current, just want t+1...t+k, so remove first value from y.future and y.pred.
	y.draw = fit$y.pred.draw
	pred.err = y.draw - y.future
	mse = sum(pred.err^2)/length(pred.err)
	
	#y.pred = fit$y.pred.mean
	#pred.err = y.pred-y.future
	#mse = sum(pred.err^2) / (length(pred.err))
	
	return(list(y.pred=y.draw,y=y.future,pred.err=pred.err,mse=mse))
}

#================================================================
# DLM Model Testing Cross-Val ===================================
#================================================================

dlm.fittest.cv = function(y,F,G,K=100,m0,C0,a.y=1,b.y=2,a.theta,b.theta,win.size=10000,iter=11000,burn=1000){
	#-------------------------------------------------------------
	#FUNCTION: 	Repeats the dlm.fittest function for a moving window of 
	#			specified size.  Calculates MSE for each window, and saves
	#			predicted vs actual values for each window's fit test.
	#-------------------------------------------------------------
	#INPUTS:	m0 = prior mean for theta.
	#			C0 = prior (diagonal) cov matrix for theta.
	#			y = vector of responses
	#			F = matrix of covariates, where each row
	#				Ft = vector of covariates at time t.
	#			F.future = matrix of covariates for times t through K.
	#			G = state matrix.
	#			a.y, b.y = Hyperparams for v ~ IG(a,b).  a=1,b=2 = IG(.5,.5).
	#			a.theta, b.theta = Vectors of hyperparams for each w_i ~ IG(a_i,b_i).
	#			iter = Number of iterations for Gibbs Sampler.
	#			burn = Number of iterations to burn for Gibbs Sampler.
	#			win.size = Number of data points in each window of 'known' observations.
	#-------------------------------------------------------------	
	#OUTPUTS:	y.pred.mean = predicted y means for t=0, t+1, t+2, ..., t+K.
	#			y.actual = actual y values for the K-step-ahead window.
	#			pred.err = y.pred.mean - y.actual for each predicted timepoint.
	#			MSE = mean squared error for range of predicted timepoints.
	#-------------------------------------------------------------
	
	#Calculate number of windows. (Windows shift forward by win.size/2 each time.)
	n = length(y)

	start = 0
	end = win.size
	if(end + K > n){
		print('Error: win.size + K cannot exceed length(y).')
	}
	
	t0 = start	#Vector to hold start points for each window.
	tn = end	#Vector to hold endpoints for each window.
	win = 1
	
	for (i in 2:n){
		start 	= t0[i] = end - win.size/2
		end 	= tn[i] = start + win.size
		if (end + K > n) break
		win = win + 1
	}
	
	#-------------------------------------------------------------
	#Initialize predicted, actual y values, error, and mse.
	y.pred = matrix(0,K,win)	#Each column is a vector of K predicted values for a single window.
	y.known = matrix(0,K,win)	#Each col is a vector of K 'known' values for a single window.
	error = matrix(0,K,win)		#Each col is a vector of K errors for a single window.
	mse = rep(0,win)			#Vector of MSE values for each window.
		
	#Loop through each window.
	for (i in 1:win){
		#Run test.
		window=c(t0[i],tn[i])
		test = dlm.fittest(y,F,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,window,iter,burn)
		
		#Save values.
		y.pred[,i] = test$y.pred
		y.known[,i] = test$y
		error[,i] = test$pred.err
		mse[i] = test$mse
	}
	
	#Return output.
	return(list(y.pred=y.pred, y.known=y.known, error=error, mse=mse))	
}