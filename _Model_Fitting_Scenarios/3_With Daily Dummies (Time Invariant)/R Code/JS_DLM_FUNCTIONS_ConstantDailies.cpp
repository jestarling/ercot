#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <cmath.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//================================================================2
// MODEL: ========================================================
//================================================================

//	y = F_{t} * theta_{t} + v, 		v ~ N(0,V)
//									
//									V ~ IG(alpha_{y},beta_{y})
//									alpha_{y} = a_{y}^2 / b_{y}
//									beta_{y} = a_{y} / b_{y}
//									s.t. a_{y} = prior mean for each w_{y}
//									and b_{y} = prior variance for each w_{y}

//	theta_{t} = G * theta_{t-1}, 	w ~ N(0,W)
//									W = diag(w1...wp)

//	w1,...,wp ~ IG(alpha_i,beta_i)
//									alpha_{i} = a_{i}^2 / b_{i}
//									beta_{i} = a_{i} / b_{i}
//									s.t. a_{i} = prior mean for each w_{i}
//									and b_{i} = prior variance for each w_{i}


//#############################################################################
//namespace dlm: support functions which do not need to be called directly in R.

namespace dlm{
	
	//================================================================
	// rcpp_seq: (replicates R's seq() fctn) =========================
	//================================================================		
	struct add_multiple {
		int incr;
		int count;
		add_multiple(int incr)
			: incr(incr), count(0)
		{}
			inline int operator()(int d) {
				return d + incr * count++;
			}
		};
	
		Rcpp::NumericVector rcpp_seq(double from_, double to_, double by_ = 1.0) {
			int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
			int from = adjust * from_;
			int to = adjust * to_;
			int by = adjust * by_;

			std::size_t n = ((to - from) / by) + 1;
			Rcpp::IntegerVector res = Rcpp::rep(from, n);
			add_multiple ftor(by);

			std::transform(res.begin(), res.end(), res.begin(), ftor);
			return Rcpp::NumericVector(res) / adjust;
		}
		
		//================================================================
		// rmvnormArma: ==================================================
		//================================================================	

		mat rmvnormArma(int n, vec mu, mat sigma) {
			//-------------------------------------------------------------
			// FUNCTION: 	Generates realizations from multivariate normal.	
			//-------------------------------------------------------------	
			int ncols = sigma.n_cols;
			mat Y = randn(n, ncols);
			mat result = (repmat(mu, 1, n).t() + Y * chol(sigma)).t();
			return result;
		}
} //end namespace dlm.



//================================================================
// Full Conditional for Time-Constant Variance V =============
//================================================================
// [[Rcpp::export]]	

double fullcond_v(double a, double b, mat F, mat theta, vec y){
			//-------------------------------------------------------------
			// FUNCTION: 	Samples from the full conditional of v|y,F,theta,a,b.
			//				(v | ...) ~ IG(a^2/b + T/2, a/b + .5 * SSy)
			//						 
			//				SSy = /sum_{t=1}^{T} (y.t - F.t %*% theta.t)^2
			//-------------------------------------------------------------
			// INPUTS:		a 		= E(1/v), prior mean for the precision.
			//				b 		= Var(1/v), prior variance for the precision.
			//				y 		= vector of y responses, y_1,...,y_t.
			//				F 		= matrix of covariate vectors; each vector is a row.
			//				theta 	= matrix of theta values; each theta_t is a col.	
			//-------------------------------------------------------------
			// OUTPUT:		v = Scalar; one sample from the full conditional of v.
			//-------------------------------------------------------------
	
			// Number of time points.
			int T = y.n_elem; 	
	
			// Calculate SSy for rate parameter.
			double SSy = 0;
			for (int t=0; t < T; ++t){
				SSy = SSy +  as_scalar(pow(y(t) - F.row(t) * theta.col(t) ,2));
			}
	
			// Convert a, b inputs to alpha, beta.
			double alpha = pow(a,2) / b;
			double beta  = a/b;
	
			// Calculate posterior shape and rate.
			double sh = .5 * T + alpha;			//gamma shape
			double rt = .5 * SSy + beta;	//gamma rate
	
			// Draw a scalar v from inverse gamma full conditional.
			// R::rgamma uses scale instead of rate; invert rate.
			double v = 1 / Rf_rgamma(sh, 1/rt);
	
			// Return output.	
			return(v); 
}



//=====================================================================
// Full Conditional for Time-Constant Diag W Elements (indep IG Priors)
//=====================================================================
// [[Rcpp::export]]	

vec fullcond_w(vec a, vec b, mat G, mat theta, vec y, vec m0){
		//-------------------------------------------------------------
		// FUNCTION: 	For W = diag(w1...wp), samples from the full conditional
		//				for each (wi | ... ) element.
		//				(wi | ...) ~ IG(a^2/b + T/2, a/b + .5 * SS.theta)
		//						
		//				SS.theta = /sum_{t=1}^{T} (theta.t.i - (G.t %*% theta.t).i )^2
		//-------------------------------------------------------------
		// INPUTS:		a 		= E(1/v), prior means vector for the precisions 1 to p.
		//				b 		= Var(1/v), prior variances vector for the precisions 1 to p.
		//				y 		= vector of y responses, y_1,...,y_t.
		//				G 		= state matrix. (Constant in time.)
		//				theta 	= matrix of theta values; each theta_t is a column.	
		//				m0 		= prior for theta: theta0 ~ N(m0,C0)
		//-------------------------------------------------------------
		// OUTPUT:		w = a vector of samples from each of the full conditionals of w1...wp.
		//-------------------------------------------------------------

		int T = y.n_elem; 		// Number of time points.
		int p = theta.n_rows;	// Number of covariates.

		// Initialize variables.
		vec alpha 		= zeros(p);
		vec beta 		= zeros(p);
		vec sh 			= zeros(p);
		vec rt 			= zeros(p);
		vec w 			= zeros(p);		// Empty vector to hold w1,...,wp.
		vec SS_theta 	= zeros(p);  

		//-------------------------------------------------------------
		// Set up SS.theta.
		//-------------------------------------------------------------

		// Create initial vector of p theta.0 values.  Set all equal to m0.
		vec theta_0 = m0;

		// SS.theta for t=1.
		SS_theta = pow(theta.col(0) - m0,2);

		// SS.theta for t=2 to T.
		for (int t=1; t<T; ++t){
			SS_theta =+ pow(theta.col(t) - theta.col(t-1), 2);
		}

		// Convert a, b inputs to alpha, beta.
		alpha = pow(a,2)/b;
		beta = a/b;

		// Calculate posterior shape and rate.
		sh = .5*T + alpha;			// Gamma shape
		rt = .5 * SS_theta + beta;	// Gamma rate

		// Draw p scalar w's from inverse gamma full conditionals.
		// R::rgamma uses scale instead of rate; invert rates.
		for(int i=0; i<p; ++i){
			w(i) = 1 / Rf_rgamma(sh(i), 1/rt(i));
		}

		// Return output.
		return(w);
}



//================================================================
// FFBS FUNCTION: ================================================
//================================================================	
// [[Rcpp::export]]	

List ffbs(vec m0, mat C0, vec y, mat F, mat G, double v, mat W){
	//-------------------------------------------------------------
	// FUNCTION: DLM FFBS (forward filter backward sampling).
	//			For retrospectively studying system: sampling theta (states) at 1:T.
	//			Treats V and W as known, time-invariant.
	//-------------------------------------------------------------
	// INPUTS:	m0 	= prior mean for theta.
	//			C0 	= prior (diagonal) cov matrix for theta.
	//			V 	= prior scalar for obs variance.
	//			W 	= State cov matrix at time t.
	//			y 	= vector of responses.
	//			Ft 	= vector of covariates at time t.
	//			Gt 	= state matrix at time t (default Im)	
	//-------------------------------------------------------------
	// OUTPUT:	theta_ff 		= theta_{1:T} matrix (cols = times) for forward filtered theta.
	//			theta_ffbs_draw = theta_{1:T} matrix (cols = times) for draws from ffbs theta.
	//			mt 				= mean vector for theta_T.
	//			Ct 				= cov matrix for theta_T.
	//			f 				= mean for y_T.
	//			Q				= var for y_T.
	//			forecast_err	= y_T - f_T	 
	//-------------------------------------------------------------	
	
	int T = y.n_elem; 	// Number of time points.
	int p = F.n_cols;	//N umber of covariates.
	
	//Set up diagonal matrix of small values to ensure psd.
	mat Sigma = eye(p,p);
	Sigma.replace(1,1e-12);
	
	//-------------------------------------------------------------
	// Initialize data structures for Kalman Filter.
	//-------------------------------------------------------------
	
	// Initialize data structures to hold mean and var for 1-step-ahead 
	// predictive of theta_{t-1} | y_{1:t-1}, ie 'prior' for theta_t.
	mat a 	= zeros(p,T);
	cube R 	= zeros(p,p,T);
	
	// Initialize data structures to hold mean and var for 1-step-ahead 
	// predictive of y_t | y_{1:t-1}.
	vec f 	= zeros(T);			// f(t) accesses scalar f at time t.
	vec Q 	= zeros(T);			// Q(t) accesses matrix Q at time t.
	
	// Initialize data structures to hold mean and var for filtering
	// distribution, theta_t | y_{1:T}.
	mat m 	= zeros(p,T);		// m(,t) accesses vector m at time t.
	cube C 	= zeros(p,p,T);		// C(,,t) accesses matrix C at time t.
	
	// Temp vec to precache R*F'*Qinv*F*R during filtering update.
	vec At = zeros(p);

	// Holds draws of theta_t from forward filter, for t=1 to T.
	mat theta_ff = zeros(p,T);
	
	//-------------------------------------------------------------
	// Initialize data structures for FFBS.
	//-------------------------------------------------------------	
	
	// Initialize data structures to hold mean and var for backward 
	// sampling mean and var at each time t.
	mat h 	= zeros(p,T);				// h(,t) accesses vector h at time t.
	cube H 	= zeros(p,p,T);				// H(,,t) accesses matrix H at time t.
	
	// Initialize data structures to hold draws of backwards sampled theta.
	mat theta_ffbs_draw = zeros(p,T);	//theta(,t) accesses vector theta at time t.
	
	// Placeholder for forecast errors e.t = y.t - f.t for each Kalman iteration.
	double et = 0;
	
	// Temp matrix to precache C.t * G.t+1' inv(R.t+1) during filtering update.
	mat CGRinv = zeros(p,p);

	//-------------------------------------------------------------
	// Forward Filter: Time t=1 Update (indexed at t=1 as 0).
	//-------------------------------------------------------------
	
	// Initialize variables and extract data at time t: y.t (scalar), and F.t (row vector of covariates at time t).
	mat Ft		= F.row(0);		// Ft is a (1xp) row vector.
	double yt 	= y(0);
	
	// Prior for theta_t.
	a.col(0) 	= G * m0;		
	R.slice(0) 	= G * C0 * G.t() + W;	
	
	// 1-step-ahead predicive of y.t | y.t:t-1.
	f(0) = as_scalar(Ft * a.col(0));
	Q(0) = as_scalar(Ft * R.slice(0) * Ft.t() + v);
	
	// Filtering dist of theta.t | y.1:t.
	et = yt - f(0);						// Precache for mt and Ct updates.
	At = R.slice(0) * Ft.t() * 1/Q(0);	// Precache for mt and Ct updates. (West/Harrison Notation Ch 4)
	
	m.col(0) = a.col(0) + At * et;
	C.slice(0) = R.slice(0) - At * Q(0) * At.t();
		
	//-------------------------------------------------------------
	// Forward Filter: Times t=2 to T.  (Indexes from 0; start at 1)
	//-------------------------------------------------------------		
	
	for (int t=1; t<T; ++t){
		
		// Extract data at time t: y.t (scalar), and F.t (row vector of covariates at time t).
		Ft	= F.row(t);
		yt 	= y(t);
		
		// Prior for theta_t.
		a.col(t) 	= G * m.col(t-1);		
		R.slice(t) 	= G * C.slice(t-1) * G.t() + W;	
		
		// 1-step-ahead predicive of y.t | y.t:t-1.
		f(t) = as_scalar(Ft * a.col(t));
		Q(t) = as_scalar(Ft * R.slice(t) * Ft.t() + v);
		
		// Filtering dist of theta.t | y.1:t.
		et = yt - f(t);						// Precache for mt and Ct updates.
		At = R.slice(t) * Ft.t() * 1/Q(t);	// Precache for mt and Ct updates. (West/Harrison Notation Ch 4)
		
		m.col(t) 	= a.col(t) + At* et;
		C.slice(t)  = R.slice(t) - At * Q(t) * At.t() + Sigma;
		
		// Draw theta_t from filtering distribution.
		theta_ff.col(t) = dlm::rmvnormArma(1,m.col(t),C.slice(t));
		
	} // End forward filter loop.
	
	//-------------------------------------------------------------
	// Backward Sample: Step 1, draw theta_T ~ N(m_T, C_T)
	//-------------------------------------------------------------	
	
	// Saving mean and var for bookkeeping, but we drew this as the last fwd filter step; can re-use.
	h.col(T-1) = m.col(T-1);
	H.slice(T-1) = C.slice(T-1);
	
	theta_ffbs_draw.col(T-1) = theta_ff.col(T-1);	
	
	//-------------------------------------------------------------
	// Backward Sample: Step 2, for t = T-1 to 0, draw 
	// theta_t ~ N(h_t, H_t) using ffbs recursion. 
	//-------------------------------------------------------------		
	for (int t=(T-2); t >= 0; --t){
		
		CGRinv 	   = C.slice(t) * G.t() * inv(R.slice(t+1));

		h.col(t)   = m.col(t) + CGRinv * (theta_ff.col(t+1) - a.col(t+1));
		H.slice(t) = C.slice(t) - CGRinv * G * C.slice(t) + Sigma;
		
		// Draw theta_t from ffbs recursion.
		theta_ffbs_draw.col(t) = dlm::rmvnormArma(1,h.col(t),H.slice(t));
	
	} // End ffbs loop.
	
	//-------------------------------------------------------------
	// Save most recent m.t and C.t values to output, as well as prediction error.
	//-------------------------------------------------------------
	vec mt = m.col(T-1);
	mat Ct = C.slice(T-1);	
	vec forecast_err = y - f;

	//-------------------------------------------------------------
	// Return function output.
	//-------------------------------------------------------------
	return Rcpp::List::create( 
		_["theta_ff"] 		 = theta_ff,
		_["theta_ffbs_draw"] = theta_ffbs_draw,
		_["mt"] 			 = mt,
		_["Ct"] 			 = Ct,
		_["m"] 				 = m,	
		_["C"] 				 = C,
		_["f"] 				 = f,
		_["Q"] 				 = Q,
		_["forecast.err"] 	 = forecast_err
	) ;		
	
} //End ffbs function.



//================================================================
// DLM K-Step-Ahead Forecast Function ============================
//================================================================
// [[Rcpp::export]]	

List forecast(vec mt, mat Ct, mat F, mat G, double v, mat W, int K){
	//------------------------------------------------------------
	// FUNCTION: 	DLM Forecasting k-steps ahead for future y values.
	//-------------------------------------------------------------
	// INPUTS:		m0 	= prior mean for theta.
	//				C0 	= prior (diagonal) cov matrix for theta.
	//				F 	= matrix of covariates at times t+1 to t+k, each row = one time.
	//				Gt 	= state matrix at time t (default Im)	
	//				k 	= forecasts t+k thetas into future, for k >= 1.
	//-------------------------------------------------------------
	// OUTPUT:		y_pred_draw	 	= Y values drawn at times t+1 to t+K.
	//				y_pred_mean		= Y predicted means at times t+1 to t+K.
	//				y_pred_var		= Y predicted vars at times t+1 to t+K.
	//				theta_pred_draw	= theta values drawn at times t+1 to t+K.
	//				theta_pred_mean	= theta predicted means at times t+1 to t+K.
	//				theta_pred_var	= theta predicted vars at times t+1 to t+K.
	//-------------------------------------------------------------
	
	int p = F.n_cols;		//Number of covariates.
	vec a0 = mt;			//Initial vector a = m0.
	mat R0 = Ct;			//Initial matrix R0 = Ct.

	//-------------------------------------------------------------
	//Initialize data structures for Forecasting Recursion.
	//-------------------------------------------------------------
	
	// Initialize data structures to hold mean and var for k-step-ahead 
	// predictive of theta_{t+k} | y_{1:T}.
	vec f = zeros(K); 			// f(t) accesses vector f at time t+k.
	vec Q = zeros(K);			// Q(t) accesses matrix Q at time t+k.
	
	// Initialize data structures to hold mean and var for 1-step-ahead 
	// predictive of y_{t+k} | y_{1:T}.
	mat a = zeros(p,K);			// a(,t) accesses vector a at time t+k.
	cube R = zeros(p,p,K);		// R(,,t) accesses matrix R at time t+k.
	
	//Initialize data structures to hold draws of theta_{t+k} and y_{t+k}.
	mat theta_tPk_draw = zeros(p,K);
	vec y_tPk_draw = zeros(K);

	//-------------------------------------------------------------
	// Recursion for k=1; time (t+1).
	//-------------------------------------------------------------	
	
	// Extract future covariate row vector.
	mat Ft = F.row(0);		//Ft is a (1xp) row vector.
		
	// Mean and var for theta_{t+k} | y_{t:T}.
	a.col(0) = G * a0;
	R.slice(0) = G * R0 * G.t() + W;
	
	// Draw theta_{t+k} | y_{t:T}.
	theta_tPk_draw.col(0) = dlm::rmvnormArma(1,a.col(0),R.slice(0));
	
	// Mean and var for y_{t+k} | y_{t:T}.
	f(0) = as_scalar(Ft * a.col(0));
	Q(0) = as_scalar(Ft * R.slice(0) * Ft.t() + v);
	
	//-------------------------------------------------------------
	// Recursion for k=2 to T;  times t+2 to t+K.
	//-------------------------------------------------------------	
	
	for (int k=1; k<K; ++k){
		// Extract future covariate row vector.
		mat Ft		= F.row(k);		// Ft is a (1xp) row vector.
		
		// Mean and var for theta_{t+k} | y_{t:T}.
		a.col(k) = G * a.col(k-1);
		R.slice(k) = G * R.slice(k-1) * G.t() + W;
	
		// Draw theta_{t+k} | y_{t:T}.
		theta_tPk_draw.col(k) = dlm::rmvnormArma(1,a.col(k),R.slice(k));
		
		// Mean and var for y_{t+k} | y_{t:T}.
		f(k) = as_scalar(Ft * a.col(k));
		Q(k) = as_scalar(Ft * R.slice(k) * Ft.t() + v);
		
	} // End forecast recursion.

	//-------------------------------------------------------------
	// Draw predicted y values for k=1 to K; outside loop for speed.
	//-------------------------------------------------------------	
	
	y_tPk_draw = randn(K) % sqrt(Q) + f;
		
	//-------------------------------------------------------------
	// Return function output.
	//-------------------------------------------------------------
	return Rcpp::List::create( 
				_["y_pred_draw"] 	 = y_tPk_draw,
				_["y_pred_mean"] 	 = f,
				_["y_pred_var"] 	 = Q,
				_["theta_pred_draw"] = theta_tPk_draw,
				_["theta_pred_mean"] = a,
				_["theta_pred_var"]  = R
	) ;	
} // End forecast function.



//================================================================
// DLM Gibbs Sampler =============================================
//================================================================
// [[Rcpp::export]]	

List gibbs(vec m0, mat C0, vec y, mat F, mat G, double a_y, double b_y,
	vec a_theta, vec b_theta, int B, int burn){
	//-------------------------------------------------------------
	// FUNCTION: 	Runs a Gibbs Sampler for all theta.t, v and W.t
	//				for times t=1 to t. (Latest available timepoint.)
	//				Includes 1-step-ahead forecast for theta.t.
	//-------------------------------------------------------------
	// INPUTS:		a 	 = E(1/v), prior mean vector for the precisions 1 to p.
	//				b 	 = Var(1/v), prior variance vector for the precisions 1 to p.
	//				y 	 = vector of y responses, y_1,...,y_t.
	//				G 	 = state matrix. (Constant in time.)
	//				m0 	 = prior for theta: theta0 ~ N(m0,C0)
	//				B 	 = number of iterations for gibbs sampler.
	//				burn = number of observations at beginning of chain to burn.
	//-------------------------------------------------------------
	// OUTPUT:		<UPDATE THIS>
	//-------------------------------------------------------------
	
	//Load required functions for gibbs sampler. (Only req for those not in DLM namespace.)
	Function ffbs("ffbs"); 
	
	//-------------------------------------------------------------	
	// Set up data structures and initialize Gibbs chain.
	//-------------------------------------------------------------	
	
	//Initialize number of time points and covariates.
	int T = y.n_elem; 				// Number of time points.
	int p = F.n_cols;				// Number of covariates.
	
	//Placeholders for storing Gibbs samplers from posteriors.
	cube theta = zeros(p,T,B);		// theta(,t,b) accesses theta_t for gibbs sample b.
	mat w 	   = zeros(p,B);		// w(,b) accesses diagonal for W at gibbs sample b.
	vec v 	   = zeros(B);			// v(b) accesses v at gibbs sample b.
	
	//Placeholders to store Gibbs Samples of mt and Ct, 
	//(The most recent ffbs mean and var for theta.)
	mat mt  = zeros(p,B);
	cube Ct = zeros(p,p,B);
	
	//Initialize gibbs sampler chain values.
	theta.slice(0).zeros();
	w.col(0).ones();
	v(0) = 1;
	
	//Vector to hold full W matrix.  Does not need to be stored.
	mat W = eye(p,p);
	
	//-------------------------------------------------------------		
	// Run Gibbs sampler.
	//-------------------------------------------------------------	
	
	for(int b=1; b<B; ++b){
		
		//Progress report. Prints every 50th iteration.
		if(b % 50 == 0){
			cout << (b+1) << endl;	
		}
		
		//Update theta using ffbs theta draws.
		W.diag() 		= w.col(b-1);	
		List foo 		= ffbs(m0,C0,y,F,G,v(b-1),W);
		theta.slice(b)	= Rcpp::as<arma::mat>(foo["theta_ffbs_draw"]);
		
		//Save mt and Ct from each theta ffbs.
		mt.col(b) 	= Rcpp::as<arma::vec>(foo["mt"]);
		Ct.slice(b) = Rcpp::as<arma::mat>(foo["Ct"]);
			
		//Update v.
		v(b) = fullcond_v(a_y,b_y,F,theta.slice(b),y);
		
		//Update w.
		w.col(b) = fullcond_w(a_theta, b_theta, G, theta.slice(b), y, m0);
		
		    //Daily dummies are in cols 28-33.  Set these to zero each time.
			w.submat(27, b, 32, b) = ones<vec>(6);
			
	}	 //end Gibbs Sampler.
	
	
	//-------------------------------------------------------------	
	// Burn beginning observations.
	//-------------------------------------------------------------	
	
	//No 'shed' equivalent for vectors; define a new burn vector for v.
	vec v_burn = zeros(B-burn);
	
	if(burn > 0){
		theta.shed_slices(0,burn-1);
		w.shed_cols(0,burn-1);
		v_burn = v.subvec(burn,B-1); 
	}  
	
	else{
		v_burn = v;
	}
	
	//-------------------------------------------------------------	
	// Calculate posterior means.
	//-------------------------------------------------------------	
	
	mat theta_pm = mean(theta,2);		//Slice means for theta.	
	mat w_pm 	 = mean(w,1); 			//Col Means for w.
	double v_pm  = mean(v_burn);
	
	//Posterior means for mt and Ct.
	vec mt_pm = mean(mt,1); 			//Col means for mt.
	mat Ct_pm = mean(Ct,2);				//Slice means for Ct.
		
	//Save last posterior mean theta_t vector to output for convenience.
	vec theta_pmt = theta_pm.col(theta_pm.n_cols-1);			
	
	//-------------------------------------------------------------
	// Return function output.
	//-------------------------------------------------------------
	return Rcpp::List::create( 
		_["theta"] 		= theta,
		_["w"] 			= w,
		_["v"] 			= v_burn,
		_["mt"] 		= mt,
		_["Ct"] 		= Ct,
		_["theta_pm"] 	= theta_pm,
		_["theta_pm_t"] = theta_pmt,
		_["w_pm"] 		= w_pm,
		_["v_pm"] 		= v_pm,
		_["mt_pm"] 		= mt_pm,
		_["Ct_pm"] 		= Ct_pm
	) ;	
				
} //End gibbs function.		