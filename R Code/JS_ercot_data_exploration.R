#ERCOT Forecasting

#Set working directory.
setwd('/Users/jennstarling/UTAustin/Research/ercot')

#Read in load data.
load = read.table('data/load_data.gz',row.names=NULL,sep=',',header=T,stringsAsFactors=F)

#Read in temp, dewpoint and windspeed data.
zone_temp = read.csv(file='data/weather_processed_by_zone/zone_temp.csv',header=T,row.names=1)
zone_dewpt = read.csv(file='data/weather_processed_by_zone/zone_dewpt.csv',header=T,row.names=1)
zone_windspd = read.csv(file='data/weather_processed_by_zone/zone_windspd.csv',header=T,row.names=1)

#Read population and business hour data.
counties = read.csv('data/county_data.csv',header=T,stringsAsFactors=F)
bushr = read.csv('data/is_bushour.csv',header=T,stringsAsFactors=F)

#-------------------------------------------------------

#---------------------------------------------------------------------
#Analysis of Missing Weather Data

#1. Match up load times and weather data times.
# (Should be none missing, due to imputation in weather processing file.)
dim(load)
dim(temp)	#Temp, dewpoint and windspeed all have same measurement times.

#Count loads missing temp times and show rows.
sum(!(load$Time %in% rownames(temp)))
load[!(load$Time %in% rownames(temp)),]

#2. Match up load times and business hour data.
sum(!(load$Time %in% bushr$datetime))