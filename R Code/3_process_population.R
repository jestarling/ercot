#ERCOT Forecasting

#Aggregates county-level population information for each year into zone-level information.

#Set working directory.
setwd('/Users/jennstarling/UTAustin/Research/ercot')

#Library for table joins.
library(sqldf)

#Read country-level population data
counties = read.csv('data/county_data.csv',header=T,stringsAsFactors=F)

#-------------------------------------------
#Aggregate population data by zone.
attach(counties)
zone_pops = aggregate(counties[,c(2,9:13)], by=list(ZONE),FUN=sum,na.rm=T)
colnames(zone_pops)[1] = 'ZONE'
colnames(zone_pops)[-1] = substr(colnames(zone_pops)[-1],4,8)

#Output populations by year for each zone.
write.csv(zone_pops,file='data/population_processed_by_zone/zone_pops.csv')

#-------------------------------------------
#Match population data by year to load times and output result.
load = read.table('data/load_data.gz',row.names=NULL,sep=',',header=T,stringsAsFactors=F)

zones = colnames(load)[-1]
zone_pops_datetime = data.frame(Time=load$Time)
zone_pops_datetime$YEAR = substr(load$Time,1,4)

coast = match(zone_pops_datetime$YEAR,colnames(zone_pops))

zone_pops[1,match(zone_pops_datetime$YEAR,colnames(zone_pops))]



coast = which(zone_pops[1,colnames(zone_pops) %in% zone_pops_datetime$YEAR])


sqldf("SELECT CustomerId, Product, State 
              FROM df1
              LEFT JOIN df2 USING(CustomerID)")

A[, c('ID', B[, 1])]

cbind(A$ID, A[names(A) %in% B$col1])


for (z in zones){
	
	zone_pops_datetime = cbind.data.frame(zone_pops_datetime,abc)
}



