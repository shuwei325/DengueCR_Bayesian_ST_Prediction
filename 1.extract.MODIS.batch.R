library(tidyverse)
library(MODISTools)
library(sf)
library(tictoc)
library(furrr)
#library(doParallel)

# Carga de localizaciones y fechas----
load('localizaciones_MODIS.RData')
# fecha_maxima <- NULL
# extrae_maximo <- function(i,producto){
#   fechas <- mt_dates(product = producto,
#                       lat = localizaciones_MODIS[i,3],
#                       lon = localizaciones_MODIS[i,4])
#   return(fechas[dim(fechas)[1],2])
# }
#fechas_Q1 <- map(1:dim(localizaciones_MODIS)[1],~extrae_maximo(.,'MOD13Q1'))
#fechas_A2 <- map(1:dim(localizaciones_MODIS)[1],~extrae_maximo(.,'MOD11A2'))

#fechas2 <- mt_dates(product = 'MOD11A2',lat = localizaciones_MODIS[12,3],lon = localizaciones_MODIS[12,4])
#fechas_1 <- c('2021-02-02','2021-05-09')
#fechas_2 <- c('2021-02-10','2021-05-17')


# Formateo de localizaciones según la función mt_batch_subset
df_MODIS <- localizaciones_MODIS %>%
#  filter(CCanton==605) %>%
  select(site_name=DTA,lat=Latitud,lon=Longitud)

n_locs <- dim(df_MODIS)[1]

# Descarga de datos:
# show('1-NDVI')
# baseNDVI <- NULL
# 
# extract_NDVI <- function(i){
#   lat_loc <- df_MODIS$lat[i]
#   lon_loc <- df_MODIS$lon[i]
#   fechas_site <- mt_dates(product = 'MOD13Q1',
#                           lat = lat_loc,
#                           lon = lon_loc)
#   
#   baseNDVI_t <- mt_subset(product = 'MOD13Q1',
#                            band = '250m_16_days_NDVI',
#                            lat = lat_loc,
#                            lon = lon_loc,
#                            start=fechas_site$calendar_date[1],
#                            end = tail(fechas_site,n=1)$calendar_date,
#                            site_name = df_MODIS$site_name[i],
#                            progress = F)
#   return(baseNDVI_t)
# }

plan(multisession, workers = 4)
#baseNDVI <- future_map_dfr(1:n_locs,.f = ~extract_NDVI(.))


show('2-EVI')
baseEVI <- NULL

extract_EVI <- function(i){
  lat_loc <- df_MODIS$lat[i]
  lon_loc <- df_MODIS$lon[i]
  fechas_site <- mt_dates(product = 'MOD13Q1',
                          lat = lat_loc,
                          lon = lon_loc)
  
  baseEVI_t <- mt_subset(product = 'MOD13Q1',
                          band = '250m_16_days_EVI',
                          lat = lat_loc,
                          lon = lon_loc,
                          start=fechas_site$calendar_date[1],
                          end = tail(fechas_site,n=1)$calendar_date,
                          site_name = df_MODIS$site_name[i],
                          progress = F)
  return(baseEVI_t)
}

baseEVI <- future_map_dfr(1:n_locs,.f = ~extract_EVI(.))


show('3-NIR')
baseNIR <- NULL

extract_NIR <- function(i){
  lat_loc <- df_MODIS$lat[i]
  lon_loc <- df_MODIS$lon[i]
  fechas_site <- mt_dates(product = 'MOD13Q1',
                          lat = lat_loc,
                          lon = lon_loc)
  
  baseNIR_t <- mt_subset(product = 'MOD13Q1',
                         band = '250m_16_days_NIR_reflectance',
                         lat = lat_loc,
                         lon = lon_loc,
                         start=fechas_site$calendar_date[1],
                         end = tail(fechas_site,n=1)$calendar_date,
                         site_name = df_MODIS$site_name[i],
                         progress = F)
  return(baseNIR_t)
}

baseNIR <- future_map_dfr(1:n_locs,.f = ~extract_NIR(.))


show('4-MIR')
baseMIR <- NULL

extract_MIR <- function(i){
  lat_loc <- df_MODIS$lat[i]
  lon_loc <- df_MODIS$lon[i]
  fechas_site <- mt_dates(product = 'MOD13Q1',
                          lat = lat_loc,
                          lon = lon_loc)
  
  baseMIR_t <- mt_subset(product = 'MOD13Q1',
                         band = '250m_16_days_MIR_reflectance',
                         lat = lat_loc,
                         lon = lon_loc,
                         start=fechas_site$calendar_date[1],
                         end = tail(fechas_site,n=1)$calendar_date,
                         site_name = df_MODIS$site_name[i],
                         progress = F)
  return(baseMIR_t)
}

baseMIR <- future_map_dfr(1:n_locs,.f = ~extract_MIR(.))



show('5-LST_Day')
baseLST_Dia <- NULL

extract_LST_Dia <- function(i){
  lat_loc <- df_MODIS$lat[i]
  lon_loc <- df_MODIS$lon[i]
  fechas_site <- mt_dates(product = 'MOD11A2',
                          lat = lat_loc,
                          lon = lon_loc)
  
  baseLST_Dia_t <- mt_subset(product = 'MOD11A2',
                         band = 'LST_Day_1km',
                         lat = lat_loc,
                         lon = lon_loc,
                         start=fechas_site$calendar_date[1],
                         end = tail(fechas_site,n=1)$calendar_date,
                         site_name = df_MODIS$site_name[i],
                         progress = F)
  return(baseLST_Dia_t)
}

baseLST_Dia <- future_map_dfr(1:n_locs,.f = ~extract_LST_Dia(.))


show('6-LST_Night')
baseLST_Noche <- NULL

extract_LST_Noche <- function(i){
  lat_loc <- df_MODIS$lat[i]
  lon_loc <- df_MODIS$lon[i]
  fechas_site <- mt_dates(product = 'MOD11A2',
                          lat = lat_loc,
                          lon = lon_loc)
  
  baseLST_Noche_t <- mt_subset(product = 'MOD11A2',
                             band = 'LST_Night_1km',
                             lat = lat_loc,
                             lon = lon_loc,
                             start=fechas_site$calendar_date[1],
                             end = tail(fechas_site,n=1)$calendar_date,
                             site_name = df_MODIS$site_name[i],
                             progress = F)
  return(baseLST_Noche_t)
}

baseLST_Noche <- future_map_dfr(1:n_locs,.f = ~extract_LST_Noche(.))

save(baseEVI,baseNIR,baseMIR,baseLST_Dia,baseLST_Noche,
     file = 'MODISData.RData')


save(baseEVI,baseNDVI,baseNIR,baseMIR,baseLST_Dia,baseLST_Noche,
     file = 'MODISData.RData')



