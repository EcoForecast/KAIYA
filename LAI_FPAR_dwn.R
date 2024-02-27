repos = "http://cran.us.r-project.org"
get.pkg <- function(pkg){
  loaded <- do.call("require",list(package=pkg))
  if(!loaded){
    print(paste("trying to install",pkg))
    install.packages(pkg,dependencies=TRUE,repos=repos)
    loaded <- do.call("require",list(package=pkg))
    if(loaded){
      print(paste(pkg,"installed and loaded"))
    } 
    else {
      stop(paste("could not install",pkg))
    }    
  }
}
# Load necessary packages
get.pkg("RCurl")
get.pkg("readr")
get.pkg("XML")
get.pkg("arrow")
get.pkg("devtools")
get.pkg("MODISTools")
get.pkg("EML")
get.pkg("cronR")
get.pkg("tidyverse")

lai_download <- function(d){
  LAI_file = paste("MODIS.LAI.",d-4,'.',d,".HARV.RData", sep='')
  file_path = file.path(getwd(), 'LAI', LAI_file)
  if(file.exists(file_path)){}
  else{
    subset <- MODISTools::mt_subset(product = "MCD15A3H",
                                    band = "Lai_500m",
                                    lat=42.5369,
                                    lon=-72.1727,
                                    start=d-4,
                                    end=d,
                                    km_lr = 1,
                                    km_ab = 1,
                                    site_name = "HARV")
    save(subset,file=file_path)
  }
}

fpar_download <- function(d){
  FPAR_file = paste("MODIS.FPAR.",d-4,'.',d,".HARV.RData", sep='')
  file_path = file.path(getwd(), 'FPAR', FPAR_file)
  if(file.exists(file_path)){}
  else{
    subset <- MODISTools::mt_subset(product = "MCD15A3H",
                                    band = "Fpar_500m",
                                    lat=42.5369,
                                    lon=-72.1727,
                                    start=d-4,
                                    end=d,
                                    km_lr = 1,
                                    km_ab = 1,
                                    site_name = "HARV")
    save(subset,file=file_path)
  }
}

lai_download(Sys.Date()-30)
fpar_download(Sys.Date()-30)
