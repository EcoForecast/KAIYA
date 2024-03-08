library(dplyr)
library(lubridate)

# function for updating the weather dataset
weather_download <- function(datetime, siteid){
  ds <- neon4cast::noaa_stage3()
  df <- ds |> filter(datetime==datetime,
                     site_id==siteid) |>
    collect()
}

# temp = weather_download("2020-08-08", "HARV")

# function for downloading the historical weather dataset
weather_historical_download <- function(startdate, siteid){
  ds <- neon4cast::noaa_stage3()
  df <- ds |> filter(datetime >= lubridate::as_datetime(startdate),
                     site_id == siteid) |>
    collect()
}

# variables included in the NOAA weather dataset
# variables <- c("surface_downwelling_longwave_flux_in_air", 
#                "surface_downwelling_shortwave_flux_in_air",
#                "precipitation_flux",                       
#                "air_pressure", 
#                "relative_humidity", 
#                "air_temperature", 
#                "northward_wind", 
#                "eastward_wind")

# temp = weather_historical_download("2024-03-01", "HARV")
