library(dplyr)
library(lubridate)

# function for updating the weather dataset
weather_download <- function(selectdate, siteid){
  ds <- neon4cast::noaa_stage3()
  df <- ds |> filter(datetime==as_datetime(selectdate),
                     site_id==siteid) |>
    collect()
  return(df)
}

# temp = weather_download("2023-08-08", "HARV")

# function for downloading the historical weather dataset
noaa_historical_download <- function(site, var, reference_date, end_date){
  ds <- neon4cast::noaa_stage3()
  # df <- ds |> filter(datetime >= lubridate::as_datetime(startdate),
  #                    site_id == siteid)
  
  historical_start_date <- as_datetime(reference_date)
  historical_end_date <- as_datetime(end_date)
  
  ds %>%
    dplyr::filter(site_id == site,
                  datetime >= historical_start_date,
                  datetime <= historical_end_date,
                  variable == var) %>%
    dplyr::select(datetime, prediction, parameter) %>%
    dplyr::group_by(datetime) %>%
    dplyr::summarize(mean_prediction = mean(prediction-273.15, na.rm = TRUE), 
                     sd_prediction = sd(prediction-273.15, na.rm = TRUE), 
                     ensemble=max(parameter)) %>%
    dplyr::select(datetime, mean_prediction, sd_prediction, ensemble) %>%
    dplyr::collect()
}

noaa_forecast_download <- function(site, var, reference_date) {
  endpoint <- "data.ecoforecast.org"
  bucket <- glue::glue("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/{reference_date}")
  s3 <- arrow::s3_bucket(bucket, endpoint_override = endpoint, anonymous = TRUE)

  ds <- arrow::open_dataset(s3)

  forecast_date <- as_datetime(reference_date)

  # Filter, select, and process data
  ds %>%
    dplyr::filter(site_id == site,
                  datetime >= forecast_date,
                  variable == var) %>%
    dplyr::select(datetime, prediction, parameter) %>%
    dplyr::group_by(datetime) %>%
    dplyr::summarize(mean_prediction = mean(prediction-273.15, na.rm = TRUE), 
                     sd_prediction = sd(prediction-273.15, na.rm = TRUE), 
                     ensemble=max(parameter)) %>%
    dplyr::select(datetime, mean_prediction, sd_prediction, ensemble) %>%
    dplyr::collect()
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
