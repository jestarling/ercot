#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <cmath.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//================================================================
// MODEL: ========================================================
//================================================================

//	y = F_{t} * theta_{t} + v, 		v ~ N(0,V)
//									where
//									V ~ IG(alpha_{y},beta_{y})
//									alpha_{y} = a_{y}^2 / b_{y}
//									beta_{y} = a_{y} / b_{y}
//									s.t. a_{y} = prior mean for each w_{y}
//									and b_{y} = prior variance for each w_{y}

//	theta_{t} = G * theta_{t-1}, 	w ~ N(0,W)
//									W = diag(w1...wp)

//	w1,...,wp ~ IG(alpha_i,beta_i), where
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
	//FUNCTION: 	Generates realizations from multivariate normal.	
	//-------------------------------------------------------------	
	   int ncols = sigma.n_cols;
	   mat Y = randn(n, ncols);
		 mat result = (repmat(mu, 1, n).t() + Y * chol(sigma)).t();
		 return result;
}
} //end namespace dlm.


//================================================================
// DLM Full Conditional for Time-Constant Variance V =============
//================================================================
// [[Rcpp::export]]	
double fullcond_v(double a, double b, mat F, mat theta, vec y){
	//-------------------------------------------------------------
	//FUNCTION: 	Samples from the full conditional of v|y,F,theta,a,b.
	//			v | ... ~ IG(a^2/b + T/2 - 1, .5 * \sum_{t=1}^{T}(y_t - F'_t theta_t)^2 + a/b)
	//-------------------------------------------------------------
	//INPUTS:	a = E(1/v), prior mean for the precision.
	//			b = Var(1/v), prior variance for the precision.
	//			y = vector of y responses, y_1,...,y_t.
	//			F = matrix of covariate vectors; each vector is a row.
	//			theta = matrix of theta values; each theta_t is a col.	
	//-------------------------------------------------------------
	//OUTPUT:	v = a sample from the full conditional of v.
	//-------------------------------------------------------------
	
	int T = y.n_elem; 	//Number of time points.
	
//	mat Ftheta = F*theta;	//temp to extract diagonal in next line.
//	int SS_y = sum(pow(y - Ftheta.diag(), 2)); //sum_{t=1}^{T}(y_t - Ft*theta_t)^2
	
	double SSy = 0;
	for (int h=0; h < T; ++h){
		SSy = SSy +  as_scalar(pow(y(h) - F.row(h) * theta.col(h) ,2));
	}
	
	double alpha = pow(a,2) / b;
	double beta = a/b;
	
	double sh = T/2 + alpha;			//gamma shape
	double rt = .5 * SSy + beta;	//gamma rate
	
	double v = 1 / R::rgamma(sh,rt);
		//as_scalar(randg(1,sh,1/rt)); //randg uses scale parameterization; invert rate.
	return(v); 
}

//==================================================================
// DLM Full Conditional for Time-Constant Diag W (Indep IG Priors) =
//==================================================================
// [[Rcpp::export]]	
vec fullcond_w(vec a, vec b, mat G, mat theta, vec y, vec m0){
	//-------------------------------------------------------------
	//FUNCTION: 	Assuming W = diag(w1...wp), samples from the full
	//			conditional of each v|y,F,theta,a,b.
	//			v | ... ~ IG(a^2/b + T/2 - 1, .5 * \sum_{t=1}^{T}(theta_t - G theta_{t-1})^2 + a/b)
	//-------------------------------------------------------------
	//INPUTS:	a = E(1/v), prior mean vector for the precisions 1 to p.
	//			b = Var(1/v), prior variance vector for the precisions 1 to p.
	//			y = vector of y responses, y_1,...,y_t.
	//			G = state matrix. (Constant in time.)
	//			theta = matrix of theta values; each theta_t is a column.	
	//			m0 = prior for theta: theta0 ~ N(m0,C0)
	//-------------------------------------------------------------
	//OUTPUT:	w = a sample from the full conditionals of w1...wp.
	//-------------------------------------------------------------
	
		int T = y.n_elem; 		//Number of time points.
		int p = theta.n_rows;	//Number of covariates.
		
		//Initialize variables.
		vec alpha = zeros(p);
		vec beta = zeros(p);
		vec sh = zeros(p);
		vec rt = zeros(p);
		vec w = zeros(p);	//empty vector to hold w1,...,wp.
		vec SS_theta = zeros(T);
		
		//SS.theta for t=1.
		SS_theta = pow(theta.col(1) - G * m0,2);	
		
		//SS.theta for t=2 to T.
		for(int t=1; t<T; ++t){
			SS_theta += pow(theta.col(t) - G * theta.col(t-1),2);
		}
		
		alpha = pow(a,2)/b;
		beta = a/b;
		
		sh = T/2 + alpha;						//gamma shape
		rt = .5 * SS_theta + beta;	//gamma rate
		
		//Generate w values.
		for(int i=0; i<p; ++i){
			w(i) = 1 / R::rgamma(sh(i),rt(i));
		}
		
		return(w);
}

//} //End namespace dlm.
//#############################################################################&=

//================================================================
// FFBS FUNCTION: ================================================
//================================================================	

// [[Rcpp::export]]	
List ffbs(vec m0, mat C0, vec y, mat F, mat G, double v, mat W, int bs){
	//-------------------------------------------------------------
	//FUNCTION: 	DLM FFBS (forward filter backward smoothing).
	//			For sampling theta_s | y_{1:t} when s < t.
	//			For retrospectively studying system: sampling states at 1:t.
	//			***Treats V and W as fixed in time and known.
	//-------------------------------------------------------------
	//INPUTS:	m0 = prior mean for theta.
	//			C0 = prior (diagonal) cov matrix for theta.
	//			V = prior scalar for obs variance.
	//			W = State cov matrix at time t.
	//			y = vector of responses.
	//			Ft = vector of covariates at time t.
	//			Gt = state matrix at time t (default Im)	
	//			bs = 0 or 1, indicates if function should filter all 
	//				thetas back in time, or just return most recent theta draw.	
	//-------------------------------------------------------------
	//OUTPUT:	m = posterior mean for theta at time t-1.
	//			C = posterior cov matrix for theta at time t-1.
	//			theta.pred.t = predicted theta value at time t.
	//			y.pred.t = predicted y at time t.
	//			pred.err = yt - y.pred.
	//-------------------------------------------------------------	
	
	int T = y.n_elem; 	//Number of time points.
	int p = F.n_cols;		//Number of covariates.
	
	//Placeholders for update variables at each time t.
	mat a 	= zeros(p,T);		//a(,t) accesses vector a at time t.
	mat m 	= zeros(p,T);		//m(,t) accesses vector m at time t.
	cube R 	= zeros(p,p,T);	//R(,,t) accesses matrix R at time t.
	cube C 	=	zeros(p,p,T);	//C(,,tt) accesses matrix C at time t.
	
	//Placeholders for backward smoothing variables at each time t.
	mat h = zeros(p,T);			//h(,t) accesses vector h at time t.
	cube B = zeros(p,p,T);	//B(,,t) accesses matrix B at time t.
	cube H = zeros(p,p,T);	//H(,,t) accesses matrix H at time t.
	
	//Placeholder for drawn theta (ff) at each time t.
	mat theta = zeros(p,T);				//theta(,t) accesses vector theta at time t.

	//-------------------------------------------------------------
	// Forward Filter: Time t=1 Update:
	mat Ft = F.row(0);
	double yt = y(0);
	a.col(0) = m0;		//Store for all time for backwards smoothing.
	R.slice(0) = C0;	//Store for all time for backwards smoothing.
	
	double ft = as_scalar(Ft * a.col(0));
	double Qt = as_scalar(v + Ft * R.slice(0) * Ft.t());
	vec At = R.slice(0) * Ft.t() / Qt;
	double et = yt - ft;
		
	//Update m and C.  These are the mean and var of theta at t=1.
	m.col(0) = a.col(0) + At * et;
	C.slice(0) = R.slice(0) - At * At.t() * Qt;	
		
	//Draw theta_0.
	theta.col(0) = dlm::rmvnormArma(1,m.col(0),C.slice(0));
	
	//-------------------------------------------------------------
	//Forward Filter: t=2 to T updates.  (Indexes from 0.)
	for(int t=1; t<T; ++t){
		 Ft = F.row(t);
		 yt = y(t);
		 a.col(t) = G * m.col(t-1);
		 R.slice(t) = G * R.slice(t-1) * G.t() + W;
		 
		 ft = as_scalar(Ft * a.col(t));
		 Qt = as_scalar(v + Ft * R.slice(t) * Ft.t());
		 At = R.slice(t) * Ft.t() / Qt;
		 et = yt - ft;
		 
		 //Update m and C.
		 m.col(t) = a.col(t) + At * et;
		 C.slice(t) = R.slice(t) - At * At.t() * Qt;	
		 
		 //Draw theta_t.
		 theta.col(t) = dlm::rmvnormArma(1,m.col(t),C.slice(t));
		 
		 //Precache values for backwards sampling. H = backwards var, h = backwards mean.
		 B.slice(t-1) = C.slice(t-1) * G.t() * inv(R.slice(t));
		 H.slice(t-1) = C.slice(t-1) - B.slice(t-1) * (R.slice(t) - B.slice(t-1)) 
				* inv(R.slice(t)) * G * C.slice(t-1);
		 h.col(t-1) = m.col(t-1) + B.slice(t-1) * (theta.col(t) - a.col(t));
		 			
	}	//end forward filter loop.
	
	//-------------------------------------------------------------
	//Backwards sampling.  (Drawing thetas, not just returning smoothed theta means.)
	
	//Variable to set column dimension of theta_draw matrix.
	int theta_cols = 1;
	if (bs ==1){
		theta_cols = T;
	}
	
	//Initialize theta_draw matrix.
	mat theta_draw = zeros(p,theta_cols);
	
	//If no backwards sampling.
	if (bs ==0){
		theta_draw = dlm::rmvnormArma(1,m.col(T-1),C.slice(T-1));
	}
	
	//If backwards sampling.
	if (bs ==1){
		theta_draw.col(T-1) = dlm::rmvnormArma(1,m.col(T-1),C.slice(T-1));
		
		for( int t = (T-2); t >= 0; --t){
			theta_draw.col(t) = dlm::rmvnormArma(1,h.col(t),H.slice(t));
		}	
	}
	
	//-------------------------------------------------------------
	//Extract most recent theta mean and variance to return.
	vec mt =  m.col(T-1);
	mat Ct = C.slice(T-1);
	
	//-------------------------------------------------------------
	//RETURN OUTPUT.
	return Rcpp::List::create( 
		_["theta_draw"] = theta_draw,
		_["mt"] = mt,
		_["Ct"] = Ct,
		_["m_all_t"] = m,	
		_["C_all_t"] = C
	) ;		
	
}	 //End ffbs function.

//================================================================
// DLM K-Step-Ahead Forecast Function ============================
//================================================================

// [[Rcpp::export]]	
List forecast(vec mt, mat Ct, mat F, mat G, double v, mat W, int K){
	//------------------------------------------------------------
	//FUNCTION: 	DLM Forecasting k-steps ahead for future y values.
	//			** If K = 0, returns Kalman Filter results at time t.
	//			** If K > 0, returns k-step-ahead forecast, including results at time t. 
	//-------------------------------------------------------------
	//INPUTS:	m0 = prior mean for theta.
	//				C0 = prior (diagonal) cov matrix for theta.
	//				V = prior scalar for obs variance.
	//				W = State cov matrix at time t.
	//				F = matrix of covariates at times t+1 to t+k, each row = one time.
	//				Gt = state matrix at time t (default Im)	
	//				k = forecasts t+k thetas into future, for k >= 1.
	//-------------------------------------------------------------
	//OUTPUT:	m = posterior mean for theta at time t-1.
	//				C = posterior cov matrix for theta at time t-1.
	//				theta.pred.t = predicted theta value at time t.
	//				y.pred.t = predicted y at times t, t+1, ..., t+K.
	//				pred.err = yt - y.pred at each time.
	//-------------------------------------------------------------
	
	int T = K; 					//Number of time points.
	int p = F.n_cols;		//Number of covariates.
	
	//Placeholders for update variables at each time t.
	mat a 	= zeros(p,T);		//a(,t) accesses vector a at time t.
	mat m 	= zeros(p,T);		//m(,t) accesses vector m at time t.
	cube R 	= zeros(p,p,T);	//R(,,t) accesses matrix R at time t.
	
	vec f = zeros(T); 			//f(,t) accesses vector f at time t.
	vec Q = zeros(T);				//Q(,,t) accesses matrix Q at time t.
	
	//--------------------------------
	// Values at time t+1.
	
	//theta mean and var updates. theta_{t+k} ~ N(a_{t+k},R_{t+k})
	a.col(0) = G * mt; 								//Because a0 = mt.
	R.slice(0) = G * Ct * G.t() + W; 	//Because R0 = Ct.
	
	// y mean and var updates.  y_{t} ~ N(ft,Qt)
	mat Ft = F.row(0);	//extract current covariates.
	f(0) = as_scalar(Ft * a.col(0));
	Q(0) = as_scalar(v + Ft * R.slice(0) * Ft.t());
	
	//--------------------------------
	//Values at times t+2, ..., t+K:
		for(int k=1; k<T; ++k){

		//theta mean and var updates. theta_{t+k} ~ N(a_{t+k},R_{t+k})
		a.col(k) = G * a.col(k-1);
		R.slice(k) = G * R.slice(k-1)*G.t() + W;

		//y mean and var updates.
		Ft = F.row(k);	//extract current covariates.
		f(k) = as_scalar(Ft * a.col(k));
		Q(k) = as_scalar(v + Ft * R.slice(k) * Ft.t());

	} //end {t+k} loop.
	
	//--------------------------------
	
	//Draw predicted y values at each time.
	// 		y_{t+k} ~ N(f_{t+k}, Q_{t+k})
	//		For speed, using r(0,1) Armadillo generator and 
	//		un-standardizing with vectors of means (f) and sdevs (sqrt(Q))
	
	vec y_pred_draw = randn(T) % sqrt(Q) + f;
	
	//Rename parameters for forecast expected value and variance.  y_{t} | y_{1,...,t-1} ~ N(f_t, Q_t).	
	vec y_pred_mean = f;				//y_pred_mean(,t) accesses vector y_pred_mean at time t.
	vec y_pred_var = 	Q; 				//y_pred_var(,,t) accesses matrix y_pred_var at time t.
		
	//--------------------------------
	//RETURN OUTPUT.
		return Rcpp::List::create( 
					_["y_pred_draw"] = y_pred_draw,
					_["y_pred_mean"] = y_pred_mean,
					_["y_pred_var"] = y_pred_var
		) ;	
}	 //End forecast function.

//================================================================
// DLM Gibbs Sampler =============================================
//================================================================

// [[Rcpp::export]]	
List gibbs(vec m0, mat C0, vec y, mat F, mat G, double a_y, double b_y,
	vec a_theta, vec b_theta, int B, int burn){
		//-------------------------------------------------------------
		//FUNCTION: 	Runs a Gibbs Sampler for all theta.t, v and W.t
		//			for times t=1 to t. (Latest available timepoint.)
		//			Includes 1-step-ahead forecast for theta.t.
		//-------------------------------------------------------------
		//INPUTS:	a = E(1/v), prior mean vector for the precisions 1 to p.
		//			b = Var(1/v), prior variance vector for the precisions 1 to p.
		//			y = vector of y responses, y_1,...,y_t.
		//			G = state matrix. (Constant in time.)
		//			m0 = prior for theta: theta0 ~ N(m0,C0)
		//			B = number of iterations for gibbs sampler.
		//			burn = number of observations at beginning of chain to burn.
		//-------------------------------------------------------------
		//OUTPUT:	w = a sample from the full conditionals of w1...wp.
		//-------------------------------------------------------------
		
		//Load required functions for gibbs sampler.
		Function ffbs("ffbs"); 
		//Functon rcpp_seq("rcpp_seq");
		//Function fullcond_w("fullcond_w");
		//Function fullcond_v("fullcond_v");
		
		int T = y.n_elem; 	//Number of time points.
		int p = F.n_cols;		//Number of covariates.
		
		//Placeholders for storing Gibbs samplers from posteriors.
		cube theta = zeros(p,T,B);	//theta(,t,b) accesses theta_t for gibbs sample b.
		mat w = zeros(p,B);					//w(,b) accesses diagonal for W at gibbs sample b.
		vec v = zeros(B);						//v(b) accesses v at gibbs sample b.
		
		//Placeholders to store Gibbs Samples of mt and Ct, 
		//the most recent ffbs mean and var for theta.
		mat mt = zeros(p,B);
		cube Ct = zeros(p,p,B);
		
		//Initialize values.
		theta.slice(0).zeros();
		w.col(0).ones();
		v(0) = 1;
		
		//Vector to hold full W matrix.  Does not need to be stored.
		mat W = eye(p,p);
		int bs = 1;	//Indicator for backwards sampling; return thetas for all t.
		
		//-------------------------------------------------------------		
		//Run Gibbs sampler.
		for(int b=1; b<B; ++b){
			
			cout << b+1 << endl; //Prints current iteration.
			
			//Update theta.
			W.diag() = w.col(b-1);	
			List foo = ffbs(m0,C0,y,F,G,v(b-1),W,bs);
			theta.slice(b)= Rcpp::as<arma::mat>(foo["theta_draw"]);
			
			//Save mt and Ct from each theta ffbs.
			mt.col(b) = Rcpp::as<arma::vec>(foo["mt"]);
			Ct.slice(b) = Rcpp::as<arma::mat>(foo["Ct"]);

			//Update v.
			v(b) = fullcond_v(a_y,b_y,F,theta.slice(b),y);
			//v(b) = dlm::fullcond_v(a_y,b_y,F,theta.slice(b),y);
			
			//Update w.
			w.col(b) = fullcond_w(a_theta, b_theta, G, theta.slice(b), y, m0);
			//w.col(b) = dlm::fullcond_w(a_theta, b_theta, G, theta.slice(b), y, m0);
		}	 //end Gibbs Sampler.
		
		//Burn beginning observations.
			theta.shed_slices(0,burn-1);
			w.shed_cols(0,burn-1);
			vec v_burn = v.subvec(burn,B-1); //Cannot resize vector v, so re-defining.
		
		//-------------------------------------------------------------	
		//Calculate posterior means.
		mat theta_pm = mean(theta,2);		//Slice means for theta.	
		mat w_pm = mean(w,1); 					//Col Means for w.
		double v_pm = mean(v_burn);
		
		//Posterior means for mt and Ct.
		vec mt_pm = mean(mt,1); 			//Col means for mt.
		mat Ct_pm = mean(Ct,2);				//Slice means for Ct.
			
		//Save last posterior mean theta_t vector to output for convenience.
		vec theta_pmt = theta_pm.col(theta_pm.n_cols-1);			
		//--------------------------------
		//RETURN OUTPUT.
		return Rcpp::List::create( 
			_["theta"] = theta,
			_["w"] = w,
			_["v"] = v_burn,
			_["mt"] = mt,
			_["Ct"] = Ct,
			_["theta_pm"] = theta_pm,
			_["theta_pm_t"] = theta_pmt,
			_["w_pm"] = w_pm,
			_["v_pm"] = v_pm,
			_["mt_pm"] = mt_pm,
			_["Ct_pm"] = Ct_pm
		) ;	
				
} //End gibbs function.		
