
## Weather forecasts
## R code to fetch forecasts from weather.gov

## [[file:~/git/research.repo/energy_datamining/load_forecast/ERCOT.org::*Weather%20forecasts][Weather\ forecasts:1]]

## Fetch 3-hour temperature and dew point forecasts from weather.gov
require(XML)

get.hourly <- function(vals, index){    #extracts temperature values from xml
    temps <- vals[index]
    field <- if(index==1) "hourly" else "dew point"
    temps <- temps[sapply(temps, function(x) any(unlist(x) == field))]
    temps <- unlist(temps[[1]][sapply(temps, names) == "value"])
    return(as.numeric(temps))
}

get.forecast <- function(lat, lon, start.time, end.time){
    url <- paste("http://graphical.weather.gov/xml/sample_products/browser_interface/ndfdXMLclient.php?lat=",
                 lat, "&lon=", lon, "&product=time-series&begin=",
                 start.time, "&end=", end.time, "&temp=temp&dew=dew", sep="")
    xml.list <- xmlToList(xmlParse(url))
    location <- as.list(xml.list[["data"]][["location"]][["point"]])
    start.valid.time <- unlist(xml.list[["data"]][["time-layout"]][names(xml.list[["data"]][["time-layout"]]) == "start-valid-time"])
    vals <- xml.list[["data"]][["parameters"]]
    temps <- get.hourly(vals, 1)
    dew.points <- get.hourly(vals, 2)
    forecast <- data.frame(as.list(location),
                          "start.valid.time" = start.valid.time,
                          "temperature" = temps,
                          "dew.point" = dew.points)
    return(forecast)
}

lat <- 30.315351                        #Austin
lon <- -97.698678                       #Austin
start.time <- "2016-11-28T00:00:00"
end.time <- "2016-11-29T20:00:00"
weather <- get.forecast(lat, lon, start.time, end.time)

## Weather\ forecasts:1 ends here
