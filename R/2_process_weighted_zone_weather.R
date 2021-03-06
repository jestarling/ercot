#Convert station-level weather data to weighted values for each zone,
#using weighting scheme and stations in current ercot methodology.

library(matrixStats)	#For weighted row means for temp, dewpoint, windspeed.

#Read in load data.  (For ordering of zone info.)
load = read.table('Data/load_info/load_data.gz',row.names=NULL,sep=',',header=T,stringsAsFactors=F)

#Read in temp, dewpoint and windspeed data.
temp = read.csv('Data/weather_processed/temperature_impute.csv',header=T,row.names=1)
dewpt = read.csv('Data/weather_processed/dewpoint_impute.csv',header=T,row.names=1)
windspeed = read.csv('Data/weather_processed/windSpeed_impute.csv',header=T,row.names=1)

#Read stations and counties data.
stations = read.table('DataRaw/station_info/station_data.gz',row.names=NULL,sep=',',header=T)
stations_key = read.csv('DataRaw/station_info/station_key_usaf_jan2017.csv',header=T)

#-------------------------------------------------------
#Construct zone temp/dewpt/windspeed data, based on weighted zone data.

#Zone weighting based on current ercot methodology.
zone_weights = data.frame(
	zone = c('NORTH','NORTH','NORTH_C','NORTH_C','NORTH_C','EAST','EAST',
	'FAR_WEST','FAR_WEST','WEST','WEST','WEST','SOUTH_C','SOUTH_C','COAST','COAST',
	'COAST','SOUTHERN','SOUTHERN','SOUTHERN'),
	
	station_id = c('KSPS','KPRX','KDFW','KACT','KMWL','KTYR',
	'KLFK','KINK','KMAF','KABI','KSJT','KJCT','KAUS','KSAT',
	'KLVJ','KGLS','KVCT','KCRP','KBRO','KLRD'),
	
	weight = c(.5,.5,.5,.25,.25,.5,.5,.5,.5,.4,.4,.2,.5,.5,.5,.3,.2,.4,.4,.2)
)

zone_weights$station_name = stations_key$STATION[match(zone_weights$station_id,stations_key$ICAO)]

#-------------------------------------------------------
#PROCESS TEMPERATURES, DEWPOINTS, & WINDSPEED BY ZONE:
#Temp is provided for each station.  Calculate wavg temp for each zone
#based on weighting in zone_weights table..

#Calculate & save avg temp for each zone at each time point.
zone_names = unique(zone_weights$zone)

#Begin with empty columns.
zone_temp = zone_dewpt = zone_windspd = data.frame(matrix("",ncol=length(zone_names),nrow=nrow(temp)))
rownames(zone_temp) = rownames(zone_dewpt) = rownames(zone_windspd) = rownames(temp)
colnames(zone_temp) = colnames(zone_dewpt) = colnames(zone_windspd) = zone_names
	
#Loop through each zone.
for (name in zone_names){
	
	#Vector of stations and respective weights in each zone.
	stns_in_zone = as.character(zone_weights$station_id[zone_weights$zone==name])
	weights = zone_weights$weight[zone_weights$zone==name]
	
	#Temps, dewpoints and windspeeds for stations at each time.
	temps = temp[,stns_in_zone]
	dewpts = dewpt[,stns_in_zone]
	winds = windspeed[,stns_in_zone]
	
	#WAVG zone temp, dewpoint and windspeed for stations in zone.
	zone_temp[,name] = rowWeightedMeans(as.matrix(temps),w=weights,na.rm=T)
	zone_dewpt[,name] = rowWeightedMeans(as.matrix(dewpts),w=weights,na.rm=T)
	zone_windspd[,name] = rowWeightedMeans(as.matrix(winds),w=weights,na.rm=T)
}

#Re-order columns to match load data.
zone.order = colnames(load)[2:9]
zone_temp = zone_temp[,zone.order]
zone_dewpt = zone_dewpt[,zone.order]
zone_windspd = zone_windspd[,zone.order]

#Output processed data.
write.csv(zone_temp,file='Data/weather_processed_by_zone/zone_temp.csv')
write.csv(zone_dewpt,file='Data/weather_processed_by_zone/zone_dewpt.csv')
write.csv(zone_windspd,file='Data/weather_processed_by_zone/zone_windspd.csv')

#------------------------------------------------------------------
#Latex table for zone weights.
zone_weights_pretty$zone = c('North','North','North Central','North Central','North Central',
	'East','East',
	'Far West','Far West','West','West','West','South Central','South Central','Coast','Coast',
	'Coast','Southern','Southern','Southern')

library(xtable)
xtable(zone_weights_pretty)