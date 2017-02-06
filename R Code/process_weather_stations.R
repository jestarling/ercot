library(lubridate)
library(reshape2)
library(readr)
library(imputeTS)
library(foreach)

####
# Read in data
####

station_list = dir('data/weather_stations/')

weather_data = foreach(i = seq_along(station_list), .combine='rbind') %do% {
	cat("Reading", station_list[i], "\n")
	path_name = paste0('data/weather_stations/', station_list[i])
	this_year = read_csv(path_name, col_types = cols(
  		Time = col_datetime(),
  		dewpoint = col_double(),
  		stationName = col_character(),
  		temperature = col_double(),
  		windSpeed = col_double()
	))
	this_year
}


####
# Melt and re-cast to get three data frames of temp, dewpoint, and windSpeed
####

all_stations_melt = melt(weather_data, id = c("Time", "stationName"),
	measure = c("temperature", "dewpoint", "windSpeed"))

# Cast as array: (days) x (Stations) x (3 variables).
# The aggregation function just takes the first observation for a
# given time/station combo, thus implicitly de-duplicating the data.
out = acast(all_stations_melt, Time ~ stationName ~ variable, fun = function(x) {x[1]})

# Split the array out as three data frames
temperature = as.data.frame(out[,,1])
dewpoint = as.data.frame(out[,,2])
windSpeed = as.data.frame(out[,,3])
#time_stamp = ymd_hms(rownames(temperature))
Time = rownames(temperature)

# Write out files
write.csv(temperature, file="data/weather_processed/temperature.csv")
write.csv(dewpoint, file="data/weather_processed/dewpoint.csv")
write.csv(windSpeed, file="data/weather_processed/windSpeed.csv")

# Remove those stations without at least 90% of the temp observations present
na_frac = apply(temperature, 2, function(x) sum(is.na(x)/length(x)))
scrub_stations = which(na_frac > 0.1)

# temperatures for the retained stations
temperature_scrub = temperature[,-scrub_stations]
dewpoint_scrub = dewpoint[,-scrub_stations]
windSpeed_scrub = windSpeed[,-scrub_stations]

# Write out scrubbed files
write.csv(temperature_scrub, file="data/weather_processed/temperature_scrubbed.csv")
write.csv(dewpoint_scrub, file="data/weather_processed/dewpoint_scrubbed.csv")
write.csv(windSpeed_scrub, file="data/weather_processed/windSpeed_scrubbed.csv")

###
# Impute the missing observations
#	(including times in load which are missing from temp, of which there are 18)
###
my_weather_smoother = function(y) {
	yrange = range(y, na.rm=TRUE)
	yhat = na.interpolation(y, option='linear')
	yhat = pmin(yrange[2], pmax(yrange[1], yhat))
}

#Read in load data.
load = read.table('data/load_data.gz',row.names=NULL,sep=',',header=T,stringsAsFactors=F)

#Check for load times which are missing weather reading times.
#(Temp, dewpoint and windspeed all have same measurement times.)
sum(!(load$Time %in% rownames(temp)))

#Add rows with new load times to imputed data frames.
new_times = load[!(load$Time %in% rownames(temp)),1]
temprows = as.data.frame(matrix(NA,nrow=length(new_times),ncol=ncol(temp)))
colnames(temprows) = colnames(temp)
rownames(temprows) = new_times

temperature_temp = rbind(temperature_scrub,temprows)
temperature_temp = temperature_temp[order(row.names(temperature_temp)),]

dewpoint_temp = rbind(dewpoint_scrub,temprows)
dewpoint_temp = dewpoint_temp[order(row.names(dewpoint_temp)),]

windSpeed_temp = rbind(windSpeed_scrub,temprows)
windSpeed_temp = windSpeed_temp[order(row.names(windSpeed_temp)),]

#Perform data imputation.
temperature_impute = apply(temperature_temp, 2, my_weather_smoother)
dewpoint_impute = apply(dewpoint_temp, 2, my_weather_smoother)
windSpeed_impute = apply(windSpeed_temp, 2, my_weather_smoother)

#Add back row names for imputed data frames.
rownames(temperature_impute) = sort(c(rownames(temperature_scrub),new_times))
rownames(dewpoint_impute) = sort(c(rownames(dewpoint_scrub),new_times))
rownames(windSpeed_impute) = sort(c(rownames(windSpeed_scrub),new_times))

write.csv(temperature_impute, file="data/weather_processed/temperature_impute.csv")
write.csv(dewpoint_impute, file="data/weather_processed/dewpoint_impute.csv")
write.csv(windSpeed_impute, file="data/weather_processed/windSpeed_impute.csv")





###
# A function for imputing NAs with an ARIMA (1,0,1) with season (1,0,1) component
### Not currently used
my_weather_smoother = function(y) {
	my_model = arima(y, order = c(1,0,1),
		seasonal = list(order=c(1,0,1), period=24))
	yhat = na.kalman(y, model=my_model$model)
	yhat
}

