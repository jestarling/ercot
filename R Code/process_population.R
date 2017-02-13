#ERCOT Forecasting

#Aggregates county-level population information for each year into zone-level information.

#Set working directory.
setwd('/Users/jennstarling/UTAustin/Research/ercot')

#Library for table joins.
library(sqldf)

#Read country-level population data
counties = read.csv('data/county_data.csv',header=T,stringsAsFactors=F)

#Aggregate population data by zone.
attach(counties)
zone_pops = aggregate(counties[,c(2,9:13)], by=list(ZONE),FUN=sum,na.rm=T)
colnames(zone_pops)[1] = 'ZONE'

#Output results:
write.csv(zone_pops,file='data/population_processed_by_zone/zone_pops.csv')

