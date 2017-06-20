#ERCOT Forecasting - Basic DLM Model

#================================================================
# Read in load and dlm data =====================================
#================================================================

#Housekeeping
rm(list=ls())

#Set working directory.
setwd('/Users/jennstarling/UTAustin/Research/ercot')

#Load functions.
source('R Code/JS_DLM_FUNCTIONS.R')
sourceCpp(file='R Code/JS_DLM_FUNCTIONS.cpp')

#Read in data.
y_all = readRDS('R Data Objects/dlm_y_allzones.rda')
F_all = readRDS('R Data Objects/dlm_F_allzones.rda')
G_all = readRDS('R Data Objects/dlm_G_allzones.rda')

n_all = unlist(lapply(y_all,nrow))	#Sample sizes for each zone.
p_all = unlist(lapply(F_all,ncol))	#Number of predictors for each zone. 

#================================================================
# Quick Test/Demo of Individual Functions =======================
#================================================================

#Pick a zone.
zone = 1	#COAST

#Extract DLM known data for selected zone.
n = n_all[[zone]]
p = p_all[[zone]]

y = y_all[[zone]][,1,drop=T] 
F = F_all[[zone]]
G = G_all[[zone]]

#De-mean and scale data.
ybar = mean(y)	#Save mean for un-scaling later.
y = c(scale(y))
F = as.matrix(cbind.data.frame(int=F[,1],scale(F[,2:4]),F[,5:p]))

#---------------------------------------------------------------
#Set up hyperparameters.
m0 = rep(0,p)
C0 = diag(1,p)

a.y = .1
b.y = .1	

a.theta = rep(.1,p)
b.theta = rep(.1,p)

K = 8	#Number of values into future to forecast.

#---------------------------------------------------------------
# Choose a subset of data, for speed.
t.known = seq(100,199,by=1)
y.known = y[t.known]
F.known = F[t.known,]
y.future = y[max(t.known+1):(max(t.known)+K)]
F.future = F[max(t.known+1):(max(t.known)+K),]

#---------------------------------------------------------------
# 1. Test Rcpp Gibbs Sampler.

test.gibbs = gibbs(m0,C0,y.known,F.known,G,a.y,b.y,a.theta, b.theta, B=10, burn=3)

test.gibbs
names(test.gibbs)

v = test.gibbs$v_pm
W = diag(as.numeric(test.gibbs$w_pm))

#TRY SETTING HOURLY DUMMIES TO 0 for W.
#W[,28:33] = 0

#---------------------------------------------------------------
# 2. Test Rcpp ffbs function.

test.ffbs = ffbs(m0,C0,y.known,F.known,G,v,W)
test.ffbs
names(test.ffbs)

mt = test.ffbs$mt
Ct = test.ffbs$Ct


#---------------------------------------------------------------
# 3. Test Rcpp forecast function.

test.f = forecast(mt,Ct,F = F.future, G=G, v=v, W=W, K)
test.f
names(test.f)

cbind.data.frame(y.future,y.draw=test.f$y_pred_draw,y.mean=test.f$y_pred_mean,y.var=test.f$y_pred_var)

#---------------------------------------------------------------
# 4. Test dlm.fit function.
test.fit = dlm.fit(y.known,F.known,F.future,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,iter=110,burn=10)
test.fit
names(test.fit)

#---------------------------------------------------------------
# 5. Test dlm.fittest function.
window=c(0,10000)
test.fittest = dlm.fittest(y,F,G,K,m0,C0,a.y=1,b.y=2,a.theta,b.theta,window,iter=5,burn=2)
test.fittest
names(test.fittest)

#---------------------------------------------------------------
# 6. Test full conditionals for v and w.

theta = test.ffbs$theta_draw
v.draw = fullcond_v(a.y,b.y,F.known,theta, y.known)
w.draw = fullcond_w(a.theta,b.theta,G,theta,y.known,m0)

cbind(colnames(F.known),w.draw)

### Full conditionals sanity checks:

#---------------------------------------
#Reproduce calc for v full conditional.
#---------------------------------------

T = ncol(theta)
alpha.y = a.y^2/b.y
beta.y = a.y / b.y

#Calculate SS.y in R.
SS.y = 0

for (t in 1:T){
	SS.y  = SS.y + (y.known[t] - F.known[t,] %*% theta[,t])^2
}

sh = alpha.y + T/2
rt = beta.y + SS.y/2

v.draw.test = 1 / rgamma(1,sh,rt)

#Compare results.
v.draw.test
v.draw

#---------------------------------------
#Reproduce calc for w full conditional.
#---------------------------------------

T = ncol(theta)
alpha.theta = a.theta^2/b.theta
beta.theta = a.theta / b.theta

#Calculate SS.theta in R.
SS.theta = (theta[,1] - m0)^2

for (t in 2:T){
	SS.theta =+ (theta[,t] - theta[,t-1])^2
}

sh = alpha.theta + T/2
rt = beta.theta + SS.theta/2
w.draw.test = 1 / rgamma(p, sh, rt)

#Compare output.
w.draw.test
as.numeric(w.draw)




















