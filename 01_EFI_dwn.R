library(neon4cast)
library(dplyr)
library(lubridate)
library(ggplot2)
library(arrow)
library(glue)
library(readr)

# function for downloading Y variable (NEE)
download_targets <- function(){
  readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_30min/terrestrial_30min-targets.csv.gz", guess_max = 1e6)
}

# function for downloading site metadata
download_site_meta <- function(){
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") 
  site_data %>% filter(as.integer(terrestrial) == 1)
}

# download target variable
targets = download_targets()
targets_nee = targets |> filter(variable=='nee')
save(targets_nee, file = 'nee.RData')

# download site metadata
site_meta = download_site_meta()
save(site_meta, file = 'site.RData')
