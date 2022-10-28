library(openxlsx)
library(readxl)
library(tidyverse)
library(dendextend)
library(ape)
library(Hmisc)
library(dtw)
library(sf)
library(NbClust)
library(dendextend)
library(gridExtra)
library(rgdal)
library(spdep)


load('./Data/datos_totales.RData')


#Archivos para el mapa
distritos_st <-st_read('./Data/distritos_shape/Distritos_de_Costa_Rica.shp') %>%
  filter(CODIGO!=60110)

names(distritos_st)
class(distritos_st)

cantones_st <- distritos_st %>%
  mutate(DTA_C = str_sub(CODIGO,start = 1,end = 3))%>%
  group_by(DTA_C)%>% summarise() %>%
  ungroup() %>% mutate(CCanton=as.numeric(DTA_C))

names(cantones_st)

canton.adj <- poly2nb(cantones_st)

W.canton <- nb2mat(canton.adj, style = "B") 
colnames(W.canton) <- cantones_st$CCanton
rownames(W.canton) <- cantones_st$CCanton

cantones <- datos_totales %>% dplyr::select(Canton,CCanton) %>%
  distinct() 
index.matrix <- match(cantones$CCanton,row.names(W.canton))
W.canton32 <- W.canton[index.matrix,index.matrix]


W.canton.rs <- nb2mat(canton.adj, style = "W") 

colnames(W.canton.rs) <- cantones_st$CCanton
rownames(W.canton.rs) <- cantones_st$CCanton

index.matrix <- match(cantones$CCanton,row.names(W.canton.rs))
W.canton.rs32 <- W.canton.rs[index.matrix,index.matrix]

save(W.canton32,W.canton.rs32,file = './Data/W.canton.vecinos.RData')
