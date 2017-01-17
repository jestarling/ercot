#ERCOT Forecasting

#Set working directory.
setwd('/Users/jennstarling/UTAustin/Research/ERCOT/Data from ARL')

#Read in data.
load = read.table('load_data.gz',row.names=NULL,sep=',',header=T)
temp_raw = read.table('temperature_all_years.gz',row.names=NULL,sep=',',header=T)

dewpt_raw = read.table('dewpoint_all_years.gz',row.names=NULL,sep=',',header=T)
stations = read.table('station_data.gz',row.names=NULL,sep=',',header=T)
counties = read.csv('county_data.csv',header=T)

#Maps each station to one of the 8 zones.
zones = read.csv('zones.csv',header=F)
colnames(zones) = c('Station','Zone')

#-------------------------------------------------------
#TEMPERATURE PROCESSING:
#Temp is provided for each station.  Calculate mean temp for each zone.

#Calculate mean zone temps at each point in time.

#Calculate & save avg temp for each zone at each time point.
zone_names = colnames(load)[2:8]

#Create empty columns in temp data frame for zone temps.
temp_raw[,zone_names] = NA

for (name in zone_names){
	
	#Vector of stations in each zone.
	stns_in_zone = as.character(zones[zones$Zone==name,1])

	#Assign row means to new column.
	temp_raw[,name] = rowMeans(temp_raw[,colnames(temp_raw) %in% stns_in_zone],na.rm=T)
}

#QUESTION:
	#NO TEMP READING FOR ERCOT.  IS ERCOT TEMP CALCULATED USING OTHERS?
	#Looks like ERCOT load is just sum of other loads.  

#Create new temp data with just zones.
temp = temp_raw[,c(1,256:262)]	#Zone temperatures.

#-------------------------------------------------------
#DEWPOINT PROCESSING:
#Dewpoint is provided for each station.  Calculate mean dewpoint for each zone.

#Create empty columns in temp data frame for zone temps.
dewpt_raw[,zone_names] = NA

for (name in zone_names){
	
	#Vector of stations in each zone.
	stns_in_zone = as.character(zones[zones$Zone==name,1])

	#Assign row means to new column.
	dewpt_raw[,name] = rowMeans(dewpt_raw[,colnames(dewpt_raw) %in% stns_in_zone],na.rm=T)
	#NO TEMP READING FOR ERCOT.  IS ERCOT TEMP CALCULATED USING OTHERS?
	#Looks like ERCOT load is just sum of other loads.  Double check this.
}

#Create new dewpt object for just the zone data.
dewpt = dewpt_raw[,c(1,256:262)]	#Zone temperatures.

#-------------------------------------------------------
#Check if Temp and Dewpoint and Load are measured at same times.
t1 = load$Time
t2 = temp$Time
t3 = dewpt$Time

length(t1)
length(t2)
length(t3)

identical(t2,t3) #Dewpoint and temp have same measurement times.

#Load measurements are at different times.
#NEXT STEP: Match up load and temp/dewpt obs!

#Update temp and dewpt objects to contain only times included in load data set.
temp = temp[temp$Time %in% load$Time,]
dewpt = dewpt[dewpt$Time %in% load$Time,]

#QUESTION:
	#temp and dewpt have 52570 rows.
	#load has 56952 rows.
	#So there are 4382 load times where temp/dewpt are missing.  How to treat these?
	
#QUESTION:
	#County data.  Has pop for a few years, and some other vars.
	#Can pull in the pop for each load measurement in time.
	#Want to do anything with other vars in here?  Talk about what they mean.	

#-------------------------------------------------------
#COUNTY INFO PROCESSING:

pop = load; pop[,2:10] = NA		#data frame with same times and col headers as load.

#Add year column to pop
pop$Year = as.numeric(substr(pop$Time,1,4))

pop_by_year = cbind(
	with(counties, aggregate(x=counties$POP2010, by=list(counties$ZONE), FUN=sum)),
	with(counties, aggregate(x=counties$POP2011, by=list(counties$ZONE), FUN=sum)),
	with(counties, aggregate(x=counties$POP2012, by=list(counties$ZONE), FUN=sum)),
	with(counties, aggregate(x=counties$POP2013, by=list(counties$ZONE), FUN=sum)),
	with(counties, aggregate(x=counties$POP2014, by=list(counties$ZONE), FUN=sum))
)

pop_by_year = pop_by_year[,c(1,2,4,6,8,10)]
colnames(pop_by_year) = c('Zone','2010','2011','2012','2013','2014')

#-------------------------------------------------------
#SAVE R DATA OBJECTS FOR EASY USE:

#Save original data.
setwd('/Users/jennstarling/UTAustin/Research/ERCOT/R Data Objects')
save(temp_raw,file='temp_raw.rda')
save(zones,file='zones.rda')
save(dewpt_raw,file='dewpt_raw.rda')
save(stations,file='stations.rda')
save(counties,file='counties.rda')

#Save cleaned temp and dewpoint data.
save(temp,file='temp.rda')
save(dewpt,file='dewpt.rda')

#Save population data.
#save(pop_by_year,file='pop_by_year.rda')

#-------------------------------------------------------
#Code to load data from R objects.
setwd('/Users/jennstarling/UTAustin/Research/ERCOT/R Data Objects')
load('temp.rda')
load('dewpt.rda')
load('zones.rda')
load('stations.rda')
load('counties.rda')
