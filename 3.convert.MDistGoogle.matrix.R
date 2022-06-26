load(file = './Data/distancia.RData')

load('./Data/localizaciones_MODIS.RData')

localizaciones_MDist <- localizaciones_MODIS %>% filter(str_ends(DTA,'01')) 

Distancias <- results$Distance
Tiempo <- results$Time

origenes <- map_chr(.x = 1:(dim(localizaciones_MDist)[1]),
                    ~paste(localizaciones_MDist$Latitud[.x],localizaciones_MDist$Longitud[.x])) 

canton.id <- tibble(CCanton=localizaciones_MDist$CCanton,local=origenes)

Distancias <- Distancias %>% left_join(canton.id, by=c("or"="local")) %>% rename(origen=CCanton) %>%
  left_join(canton.id, by=c("de"="local")) %>% rename(destino=CCanton)

W.Distancias <- Distancias %>% select(Distance, origen, destino) %>%
  spread(destino, Distance)
row.names(W.Distancias)<- W.Distancias$origen
W.Distancias <- as.matrix(W.Distancias[,-1])

Tiempo <- Tiempo %>% left_join(canton.id, by=c("or"="local")) %>% rename(origen=CCanton) %>%
  left_join(canton.id, by=c("de"="local")) %>% rename(destino=CCanton)

W.Tiempo <- Tiempo %>% select(Time, origen, destino) %>%
  spread(destino, Time)
row.names(W.Tiempo)<- W.Tiempo$origen
W.Tiempo <- as.matrix(W.Tiempo[,-1])

save(results,file = './Data/distancia.RData')



