library(gmapsdistance)
library(tidyverse)


load('./Data/localizaciones_MODIS.RData')

localizaciones_MDist <- localizaciones_MODIS %>% filter(str_ends(DTA,'01')) 

set.api.key(key = '')

origenes <- map_chr(.x = 1:(dim(localizaciones_MDist)[1]),
        ~paste(localizaciones_MDist$Latitud[.x],localizaciones_MDist$Longitud[.x])) 

destinos <- origenes

results <- gmapsdistance(origin = origenes,destination = destinos,
                         combinations = 'all',mode = 'driving',shape = 'long')

save(results,file = './Data/distancia.RData')
