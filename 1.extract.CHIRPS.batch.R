library(tidyverse)
library(chirps)
library(sf)
library(tictoc)

# Carga de localizaciones y fechas----

load('./Data/localizaciones_MODIS.RData')
#fechas1 <- mt_dates(product = 'MOD13Q1',lat = localizaciones_MODIS[12,3],lon = localizaciones_MODIS[12,4])
#fechas2 <- mt_dates(product = 'MOD11A2',lat = localizaciones_MODIS[12,3],lon = localizaciones_MODIS[12,4])
fechas_1 <- c('2000-02-18','2021-01-17')
fechas_2 <- c('2000-02-18','2021-02-02')


base_lat_lon <- localizaciones_MODIS %>% select(DTA,Latitud,Longitud) %>% 
  st_as_sf(coords = c('Longitud','Latitud'),crs=4326) %>%
  mutate(DTA_st = as.character(DTA))

datos_precipitacion <- get_chirps(base_lat_lon, dates = fechas_1)

save(datos_precipitacion,file = 'Data/datos_precipitacion.RData')