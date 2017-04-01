### ERCOT Forecasting - Basic DLM Model
### FUNCTIONS
### Jennifer Starling
### March 2017

library(Rcpp)
library(RcppArmadillo)


#Load C++ file.
#sourceCpp(file='/Users/jennstarling/UTAustin/Research/ercot/R Code/JS_ercot_basic_dlm_FUNCTIONS.cpp')


#================================================================
# MODEL: ========================================================
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
# DLM Kalman Filter (1-step) Function ===========================
#================================================================

dlm.Filter = function(m0,C0,Vt,W,yt,Ft,Gt){
	#-------------------------------------------------------------
	#FUNCTION: 	DLM Forward Filter based on Kalman filter.
	#			Predicts current parameter theta_t based on y_t.
	#			Also provides 1-step-ahead forecast for y_t.
	#-------------------------------------------------------------
	#INPUTS:	m0 = prior mean for theta.
	#			C0 = prior (diagonal) cov matrix for theta.
	#			V = prior scalar for obs variance.
	#			W = State cov matrix at time t.
	#			yt = scalar response at time t.
	#			Ft = vector of covariates at time t.
	#			Gt = state matrix at time t (default Im)		
	#-------------------------------------------------------------
	#OUTPUT:	m = posterior mean for theta at time t-1.
	#			C = posterior cov matrix for theta at time t-1.
	#			theta.pred.t = predicted theta value at time t.
	#			y.pred.t = predicted y at time t.
	#			pred.err = yt - y.pred.
	#-------------------------------------------------------------
	require(mvtnorm)
	
	n = length(yt)	#Dimension of y.  Should be 1 here.
	p = length(Ft)	#Number of covariates.
	
	#Update variables.
	a.t = Gt %*% m0
	f.t = crossprod(Ft,a.t)
	e.t = yt - f.t
	
	#R = G %*% C.old %*% t(G) + W	#Simplifying since Gt = G = I for all t.
	R.t = C0 + W
	Q.t = as.numeric(1 + t(Ft) %*% R.t %*% Ft)
	A.t = R.t %*% Ft / Q.t
	
	#Update m and C.
	m.t = a.t + A.t %*% e.t
	C.t = R.t - A.t %*% t(A.t) * Q.t
	
	#Store predicted y value and prediction error/var.
	y.pred.t = f.t
	y.pred.var = Q.t
	pred.err.t = e.t
	
	#Calculate predicted theta value at time t.
	theta.draw.t = rmvnorm(1,m.t,C.t) 
		
	return(list(m=m.t,
		C=C.t,
		theta.draw.t = theta.draw.t, 
		y.pred.t = y.pred.t,
		pred.err = e.t
		))
} #end dlm.Filter function

#================================================================
# DLM K-Step-Ahead Forecast Function ============================
#================================================================

dlm.kForecast = function(m0,C0,yt,V,W,F,G,K=1){
	#------------------------------------------------------------
	#FUNCTION: 	DLM Forecasting k-steps ahead.
	#			Provide forecast for y_{t+k} | y_{1:t},
	#			ie y_s | y_t when s>t.
	#-------------------------------------------------------------
	#INPUTS:	m0 = prior mean for theta.
	#			C0 = prior (diagonal) cov matrix for theta.
	#			V = prior scalar for obs variance.
	#			W = State cov matrix at time t.
	#			F = matrix of covariates at times t+1 to t+k, each row = one time.
	#			Gt = state matrix at time t (default Im)	
	#			k = forecasts t+k thetas into future, for k >= 1.
	#-------------------------------------------------------------
	#OUTPUT:	m = posterior mean for theta at time t-1.
	#			C = posterior cov matrix for theta at time t-1.
	#			theta.pred.t = predicted theta value at time t.
	#			y.pred.t = predicted y at time t.
	#			pred.err = yt - y.pred.
	#-------------------------------------------------------------
	
	require(mvtnorm)
	p = ncol(F)	#Number of covariates.
	
	#Placeholders for update variables at each time.
	a = matrix(0,nrow=p,ncol=K)	#a[,t] accesses vector a at time t.
	R = array(0,dim=c(p,p,K))	#R[,,t] accesses R matrix at time t.
	
	f = rep(0,K)	#m[,t] accesses vector m at time t.
	Q = rep(0,K)	#Q[p] access variance of y at time t.
	
	theta.pred.mean = matrix(0,nrow=p,ncol=K)
	theta.pred.var = array(0,dim=c(p,p,K))
	theta.draw = matrix(0,nrow=p,ncol=K)
	y.pred = rep(0,K)

	#------------------------------
	#INITIAL VALUES AT TIME t, ie k=0.
	
	#Distribution of theta_{t+k} ~ N(a_t(k),R_t(k)).
	a[,1] = G %*% m0
	R[,,1] = G %*% C0 %*% t(G) + W	
	
	#Distribution of y_{t+k} ~ N(f_t(k),Q_t(k)). NOTE: F_{t+k} being treated as known.
	f[1] = F[1,] %*% a[,1]
	Q[1] = t(F[1,]) %*% R[,,1] %*% F[1,] + V
	
	theta.pred.mean[,1] = a[,1]
	theta.pred.var[,,1] = R[,,1]
	y.pred[1] = F[1,] %*% theta.draw[,1] + v
	
	#------------------------------
	#Loop through k=2 to K.
	for (k in 2:K){
		#Distribution of theta_{t+k} ~ N(a_t(k),R_t(k)).
		a[,k] = G %*% a[,k-1]
		R[,,k] = G %*% R[,,k-1] %*% t(G) + W	
	
		#Distribution of y_{t+k} ~ N(f_t(k),Q_t(k)). NOTE: F_{t+k} being treated as known.
		f[k] = F[k,] %*% a[,k]
		Q[k] = t(F[k,]) %*% R[,,k] %*% F[k,] + V
		
		#Predicted values and theta draw for forecast theta_{t+k} values.
		theta.pred.mean[,k] = a[,k]	
		theta.pred.var[,,k] = R[,,k]
		theta.draw[,k] = rmvnorm(1,theta.pred.mean[,k],theta.pred.var[,,k])
		
		y.pred[k] = F[k,] %*% theta.draw[,k] + v
	}
	
	colnames(theta.pred.mean) = paste("t+",1:K,sep="")
	
	#Save mean and var for forecast y_{t+k} values.	
	y.pred.mean = f
	y.pred.var = Q
	
	return(list(theta.pred.mean=theta.pred.mean,
		theta.pred.var=theta.pred.var,
		theta.draw = theta.draw,
		y.pred = y.pred,
		y.pred.mean=y.pred.mean,
		y.pred.var=y.pred.var
		))
} #end dlm.kForecast function

#================================================================
# DLM FFBS Function =============================================
#================================================================

dlm.FFBS = function(m0,C0,y,F,G,V,W,bs=0){
	#-------------------------------------------------------------
	#FUNCTION: 	DLM FFBS (forward filter backward smoothing).
	#			For sampling theta_s | y_{1:t} when s < t.
	#			For retrospectively studying system: sampling states at 1:t.
	#			***Treats V and W as fixed in time and known.
	#-------------------------------------------------------------
	#INPUTS:	m0 = prior mean for theta.
	#			C0 = prior (diagonal) cov matrix for theta.
	#			V = prior scalar for obs variance.
	#			W = State cov matrix at time t.
	#			y = vector of responses.
	#			Ft = vector of covariates at time t.
	#			Gt = state matrix at time t (default Im)	
	#			bs = 0 or 1, indicates if function should filter all 
	#				thetas back in time, or just return most recent theta draw.	
	#-------------------------------------------------------------
	#OUTPUT:	m = posterior mean for theta at time t-1.
	#			C = posterior cov matrix for theta at time t-1.
	#			theta.pred.t = predicted theta value at time t.
	#			y.pred.t = predicted y at time t.
	#			pred.err = yt - y.pred.
	#-------------------------------------------------------------
	
	require(mvtnorm)
	
	T = length(y)	#number of time-points.
	print(T)
	p = ncol(F)		#Number of predictors.
	n = 1 			#Number of dimensions. (length y at each time t)
	
	#Placeholders for update variables at each time.
	a = matrix(0,nrow=p,ncol=T)	#a[,t] accesses vector a at time t.
	R = array(0,dim=c(p,p,T))	#R[,,t] accesses R matrix at time t.
	m = matrix(0,nrow=p,ncol=T)	#m[,t] accesses vector m at time t.
	C = array(0,dim=c(p,p,T))	#C[,,t] accesses C matrix at time t.
	
	#Placeholders for backward smoothing variables at each time.
	h =  matrix(0,nrow=p,ncol=T)	#h[,t] accesses vector h at time t.
	B = array(0,dim=c(p,p,T))		#B[,,t] accesses H matrix at time t.
	H = array(0,dim=c(p,p,T))		#H[,,t] accesses H matrix at time t.
	
	#Placeholder to hold theta.t values.
	theta = matrix(0,nrow=p,ncol=T)	#theta[,t] accesses vector theta at time t.
	
	#Placeholder for predicted y values and variances.
	y.pred = rep(0,T)
	y.pred.var = rep(0,T)
	
	#-------------------------------------------------------------
	#Forward Filter: TIME t=1 UPDATE:
	Ft = F[1,]
	yt = y[1]
	a[,1] = m0		#Store for backward smoothing.
	R[,,1] = C0		#Store for backward smoothing.
	
	f.t = crossprod(Ft,a[,1])
	Q.t = as.numeric(1 + t(Ft) %*% R[,,1] %*% Ft)
	A.t = R[,,1] %*% Ft / Q.t
	e.t = yt - f.t

	#Update m and C.
	m[,1] = a[,1] + A.t %*% e.t
	C[,,1] = R[,,1] - A.t %*% t(A.t) * Q.t
	
	#Draw a theta value.
	theta[,1] = rmvnorm(1,m[,1],C[,,1])

	#-------------------------------------------------------------
	#Forward Filter: TIMES t=2:T UPDATES:
	for (t in 2:T){
		print(t)
		Ft = F[t,]
		yt = y[t]
		a[,t] = G %*% a[,t-1]					#Store for backward smoothing.
		R[,,t] = G %*% R[,,t-1] %*% t(G) + W		#Store for backward smoothing.
	
		f.t = crossprod(Ft,a[,t])
		Q.t = as.numeric(1 + t(Ft) %*% R[,,t] %*% Ft)
		A.t = R[,,t] %*% Ft / Q.t
		e.t = yt - f.t

		#Update m and C.
		m[,t] = a[,t] + A.t %*% e.t
		C[,,t] = R[,,t] - A.t %*% t(A.t) * Q.t
		
		#Draw theta.
		theta[,t] = rmvnorm(1,m[,t],C[,,t])
		
		#Calculations for backwards smoother.
		#Update B and H. B is just precaching CG'Rinv
		B[,,t-1] = C[,,t-1] %*% t(G) %*% solve(R[,,t])
		H[,,t-1] = C[,,t-1]-B[,,t-1] %*% G %*% C[,,t-1]			#Backward var.
		h[,t-1] = m[,t-1] + B[,,t-1] %*% (theta[,t] - a[,t]) 	#Backward mean.	
    
	} #end forward filter loop.

	#-------------------------------------------------------------
	#Backward Sampling: ITERATE from T-1 to 0.
    
    #If not backwards sampling, draw only theta at time t.
    if (bs==0){
    	theta.draw = rmvnorm(1,m[,T],C[,,T])
    }
    
    #If performing backwards sampling, proceed.
    if(bs==1){
    	
    	#Set up smoothed theta matrix and first theta value.
    	theta.draw = matrix(0,nrow=p,ncol=T)	#theta.bs[,t] accesses vector theta.bs at time t-1.
    	theta.draw[,T] = rmvnorm(1,m[,T],C[,,T])
    
    	#Backwards smoothing.
    	for (t in (T-1):1){
    		print(t)
    		theta.draw[,t] = rmvnorm(1,h[,t],H[,,t])
    	}
    } #end backwards smoothing
    
  	return(list(theta.draw=theta.draw))
}


#================================================================
# DLM Full Conditional for Time-Constant Variance V =============
#================================================================

fullcond.v = function(a,b,F,theta,y){
	#-------------------------------------------------------------
	#FUNCTION: 	Samples from the full conditional of v|y,F,theta,a,b.
	#			v | ... ~ IG(a^2/b + T/2 - 1, .5 * \sum_{t=1}^{T}(y_t - F'_t theta_t)^2 + a/b)
	#-------------------------------------------------------------
	#INPUTS:	a = E(1/v), prior mean for the precision.
	#			b = Var(1/v), prior variance for the precision.
	#			y = vector of y responses, y_1,...,y_t.
	#			F = matrix of covariate vectors; each vector is a row.
	#			theta = matrix of theta values; each theta_t is a row.	
	#-------------------------------------------------------------
	#OUTPUT:	v = a sample from the full conditional of v.
	#-------------------------------------------------------------
	T = length(y)
	
	SS.y = 0
	for (t in 1:T){
		SS.y = SS.y + (y - crossprod(F[t,],theta[,t]))^2
	}
	
	#SS.y = sum((y - diag(crossprod(t(F),theta)))^2)
	
	alpha = a^2/b
	beta = a/b
	
	sh = T/2 + alpha
	rt = .5 * SS.y + beta
	
	v = 1 / rgamma(1,shape=sh,rate=rt)
	return(v)
}

#==================================================================
# DLM Full Conditional for Time-Constant Diag W (Indep IG Priors) =
#==================================================================

fullcond.w = function(a,b,G,theta,y,m0){
	#-------------------------------------------------------------
	#FUNCTION: 	Assuming W = diag(w1...wp), samples from the full
	#			conditional of each v|y,F,theta,a,b.
	#			v | ... ~ IG(a^2/b + T/2 - 1, .5 * \sum_{t=1}^{T}(theta_t - G theta_{t-1})^2 + a/b)
	#-------------------------------------------------------------
	#INPUTS:	a = E(1/v), prior mean vector for the precisions 1 to p.
	#			b = Var(1/v), prior variance vector for the precisions 1 to p.
	#			y = vector of y responses, y_1,...,y_t.
	#			G = state matrix. (Constant in time.)
	#			theta = matrix of theta values; each theta_t is a column.	
	#			m0 = prior for theta: theta0 ~ N(m0,C0)
	#-------------------------------------------------------------
	#OUTPUT:	w = a sample from the full conditionals of w1...wp.
	#-------------------------------------------------------------
	T = length(y)
	p = nrow(theta)
	
	SS.theta = rep(0,nrow(theta))
	
	#SS.theta for t=1.
	SS.theta = (theta[,1] - G %*% m0)^2
	
	#SS.theta for t=2 to T.
	for (t in 2:T){
		SS.theta = SS.theta + (theta[,t] - G %*% theta[,t-1])^2
	}
	
	alpha = a^2/b
	beta = a/b
	
	sh = T/2 + alpha
	rt = .5 * SS.theta + beta
	
	w = 1 / rgamma(p,shape=sh,rate=rt)
	return(w)
}

#================================================================
# DLM Gibbs Sampler =============================================
#================================================================

dlm.Gibbs = function(m0,C0,y,F,G,a.y,b.y,a.theta,b.theta,B,thin=0,burn=0){
	#-------------------------------------------------------------
	#FUNCTION: 	Runs a Gibbs Sampler for all theta.t, v and W.t
	#			for times t=1 to t. (Latest available timepoint.)
	#			Includes 1-step-ahead forecast for theta.t.
	#-------------------------------------------------------------
	#INPUTS:	a = E(1/v), prior mean vector for the precisions 1 to p.
	#			b = Var(1/v), prior variance vector for the precisions 1 to p.
	#			y = vector of y responses, y_1,...,y_t.
	#			G = state matrix. (Constant in time.)
	#			m0 = prior for theta: theta0 ~ N(m0,C0)
	#			B = number of iterations for gibbs sampler.
	#			thin = number of observations to thin out.
	#			burn = number of observations at beginning of chain to burn.
	#-------------------------------------------------------------
	#OUTPUT:	w = a sample from the full conditionals of w1...wp.
	#-------------------------------------------------------------
	p = ncol(F)		#Number of parameters.
	T = length(y)	#Total number of timepoints.
	
	#Array of t vectors of theta at each iteration of Gibbs sampler.
	theta = array(0,dim=c(p,T,B))  #At Gibbs iter B, each theta_t time point is a col.
	
	#Matrix to store diagonals for w at each iteration of Gibbs sampler.
	w = matrix(0,nrow=p,ncol=B)
	
	#Vector to store v at each iteration of gibbs sampler.
	v = rep(0,B)
	
	#Initialize theta.
	theta[,,1] = matrix(1,nrow=p,ncol=T,byrow=F)	#This is a bit silly - doing all theta_t initialized to m0.	
													#Doesn't matter.  Gibbs.
	#Initialize v and w.
	w[,1] = 1/rgamma(p,a.theta^2/b.theta, a.theta/b.theta)
	v[1] = 1/rgamma(1,a.y^2/b.y, a.y/b.y)
	
	#Gibbs sampler.
	for (b in 2:B){
		
		#Update theta.
		theta[,,b] = dlm.FFBS(m0,C0,y,F,G,v[b],W=diag(w[,b]))$theta.draw
		
		#Update v.
		v[b] = fullcond.v(a.y,b.y,F,theta[,,b],y)
		
		#Update w.
		w[,b] = fullcond.w(a.theta,b.theta,G,theta[,,b],y,m0)
	}
	
	#Burn beginning observations.
	if (burn > 0){
		theta = theta[,,-(1:burn)]
		v = v[-(1:burn)]
		w = w[,-(1:burn)]
	}
	
	#Thin remaining observations.
	if (thin > 0){
		theta = theta[,,seq(1,dim(theta)[3],by=thin)]
		v = v[seq(1,length(v),by=thin)]
		w = w[,seq(1,ncol(w),by=thin)]
	}
	
	#Calculate posterior means.
	theta.post.mean = apply(theta,c(1,2),mean)
	theta.post.var = apply(theta,c(1,2),var)
	
	v.post.mean = mean(v)
	v.post.var = var(v)
	
	w.post.mean = apply(w,1,mean)
	w.post.var = apply(w,1,var)
	
	#Output sampler result.
	return(list(theta=theta,v=v,w=w,
		theta.post.mean=theta.post.mean,
		theta.post.var=theta.post.var,
		v.post.mean=v.post.mean,
		v.post.var=v.post.var,
		w.post.mean=w.post.mean,
		w.post.var=w.post.var
		))
}

#================================================================
# DLM Gibbs Sampler Using CPP Functions =========================
#================================================================

dlm.gibbs = function(m0,C0,y,F,G,a.y,b.y,a.theta,b.theta,B,thin=0,burn=0){
	#-------------------------------------------------------------
	#FUNCTION: 	Runs a Gibbs Sampler for all theta.t, v and W.t
	#			for times t=1 to t. (Latest available timepoint.)
	#			Includes 1-step-ahead forecast for theta.t.
	#-------------------------------------------------------------
	#INPUTS:	a = E(1/v), prior mean vector for the precisions 1 to p.
	#			b = Var(1/v), prior variance vector for the precisions 1 to p.
	#			y = vector of y responses, y_1,...,y_t.
	#			G = state matrix. (Constant in time.)
	#			m0 = prior for theta: theta0 ~ N(m0,C0)
	#			B = number of iterations for gibbs sampler.
	#			thin = number of observations to thin out.
	#			burn = number of observations at beginning of chain to burn.
	#-------------------------------------------------------------
	#OUTPUT:	w = a sample from the full conditionals of w1...wp.
	#-------------------------------------------------------------
	p = ncol(F)		#Number of parameters.
	T = length(y)	#Total number of timepoints.
	
	#Array of t vectors of theta at each iteration of Gibbs sampler.
	theta = array(0,dim=c(p,T,B))  #At Gibbs iter B, each theta_t time point is a col.
	
	#Matrix to store diagonals for w at each iteration of Gibbs sampler.
	w = matrix(1,nrow=p,ncol=B)
	
	#Vector to store v at each iteration of gibbs sampler.
	v = rep(1,B)
	
	#Initialize theta.
	theta[,,1] = matrix(1,nrow=p,ncol=T,byrow=F)	#This is a bit silly - doing all theta_t initialized to m0.	
													#Doesn't matter.  Gibbs.
	#Initialize v and w.
	w[,1] = 1/rgamma(p,a.theta^2/b.theta, a.theta/b.theta)
	v[1] = 1/rgamma(1,a.y^2/b.y, a.y/b.y)
	
	#Gibbs sampler.
	for (b in 2:B){
		
		print(b);
		
		#Update theta.
		theta[,,b] = ffbs(m0,C0,y,F,G,v[b],W=diag(w[,b]),bs=1)$theta_draw
		
		#Update v.
		v[b] = fullcond.v(a.y,b.y,F,theta[,,b],y)
		
		#Update w.
		w[,b] = fullcond.w(a.theta,b.theta,G,theta[,,b],y,m0)
	}
	
	#Burn beginning observations.
	if (burn > 0){
		theta = theta[,,-(1:burn)]
		v = v[-(1:burn)]
		w = w[,-(1:burn)]
	}
	
	#Thin remaining observations.
	if (thin > 0){
		theta = theta[,,seq(1,dim(theta)[3],by=thin)]
		v = v[seq(1,length(v),by=thin)]
		w = w[,seq(1,ncol(w),by=thin)]
	}
	
	#Calculate posterior means.
	theta.post.mean = apply(theta,c(1,2),mean)
	theta.post.var = apply(theta,c(1,2),var)
	
	v.post.mean = mean(v)
	v.post.var = var(v)
	
	w.post.mean = apply(w,1,mean)
	w.post.var = apply(w,1,var)
	
	#Output sampler result.
	return(list(theta=theta,v=v,w=w,
		theta.post.mean=theta.post.mean,
		theta.post.var=theta.post.var,
		v.post.mean=v.post.mean,
		v.post.var=v.post.var,
		w.post.mean=w.post.mean,
		w.post.var=w.post.var
		))
}

