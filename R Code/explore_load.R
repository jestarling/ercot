library(lubridate)
library(splines)
library(imputeTS)
library(mosaic)
library(foreach)


# Grid load and weather data
load_data = read.csv("../data/load_data.csv")
temperature_impute = read.csv("../data/weather_processed/temperature_impute.csv", row.names=1)
dewpoint_impute = read.csv("../data/weather_processed/temperature_impute.csv", row.names=1)
windSpeed_impute = read.csv("../data/weather_processed/windSpeed_impute.csv", row.names=1)

# Keep the load data for dates when we have weather data
mysub = which(ymd_hms(load_data$Time) %in% ymd_hms(rownames(temperature_impute)))
load_data = load_data[mysub,]

# De-dupe the weather data by merging on first match
temp_ind = match(ymd_hms(load_data$Time), ymd_hms(rownames(temperature_impute)))
temperature_impute = temperature_impute[temp_ind,]
dewpoint_impute = dewpoint_impute[temp_ind,]
windSpeed_impute = windSpeed_impute[temp_ind,]


time_stamp = ymd_hms(load_data$Time)
all(time_stamp ==  ymd_hms(rownames(temperature_impute)))
all(time_stamp ==  ymd_hms(rownames(dewpoint_impute)))
all(time_stamp ==  ymd_hms(rownames(windSpeed_impute)))



# Weather principal components
X_all = cbind(temperature_impute, dewpoint_impute, windSpeed_impute)
pc_weather = princomp(X_all)


dim(load_data)

summary(temp_data)


plot(time_stamp, temp_interp[,1], main='Temperature')

plot(tail(time_stamp, -1), diff(temp_interp[,150])/as.numeric(diff(time_stamp)), main='Temperature fluctuations')
identify(tail(time_stamp, -1), diff(temp_interp[,150]), n=10)


my_samp = 9400:9700
plot(time_stamp[my_samp], temp_interp[my_samp,227], main='Temperature, Travis County')


pc1 = prcomp(temp_interp)
n_comps = 10
pctemp2 = pc1$x[,1:n_comps]

X_temp = foreach(i=1:n_comps, .combine='cbind') %do% {
		bs(pctemp2[,i], df=6, degree=3)
	}				
X_hour = factor(hour(time_stamp))
X_day = factor(wday(time_stamp))
X_trend = 1:nrow(load_data)


plot(pctemp2[,1], load_data$ERCOT)

lm1 = lm(load_data$ERCOT ~ X_temp + X_hour + X_day + X_trend + X_hour:X_temp)
summary(lm1)
anova(lm1)
plot(fitted(lm1), type='l')

plot(fitted(lm1), load_data$ERCOT)

## Log scale
plot(log(load_data$ERCOT))

lm2 = lm(log(load_data$ERCOT) ~ X_temp + X_hour + X_day + X_trend + X_hour:X_temp)
summary(lm2)

cor(fitted(lm1), load_data$ERCOT)
cor(exp(fitted(lm2)), load_data$ERCOT)

plot(resid(lm1), type='l')
acf(resid(lm1), lag.max=24*14); abline(v=c(24, 24*7))


# Quick and easy AR model
eps = resid(lm1)
keep_lags = c(1:4, 11:13, 22:26, 46:50, 72, 165:171)
max_lag = max(keep_lags)
start_ind = max_lag + 1
end_ind = length(eps)
train_ind = start_ind:end_ind
y = eps[start_ind:end_ind]
X_eps = foreach(lag = keep_lags, .combine='cbind') %do% {
	eps[start_ind:end_ind - lag]
}
colnames(X_eps) = paste0("lag", keep_lags)

lm_eps = lm(y ~ X_eps)
summary(lm_eps)

acf(resid(lm_eps), lag.max=24*14)

load_scrub = load_data$ERCOT[train_ind]
load_hat = fitted(lm1)[train_ind] + fitted(lm_eps)
plot(load_hat, load_scrub)
outliers = identify(load_hat, load_scrub, n=5)

test_ind = c()
for(i in outliers) {
	test_ind = c(test_ind, i + (-10:10))
}
test_ind = sort(unique(test_ind))
cbind(fitted(lm1)[train_ind][test_ind], load_hat[test_ind], load_data$ERCOT[test_ind], (load_hat - load_scrub)[test_ind])

