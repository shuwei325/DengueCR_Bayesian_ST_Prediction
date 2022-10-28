library(tidyverse)
library(sf)
#Archivos para el mapa
distritos_st <-st_read('./Data/distritos_shape/Distritos_de_Costa_Rica.shp') %>%
  filter(CODIGO!=60110)

cantones_st <- distritos_st %>%
  mutate(DTA_C = str_sub(CODIGO,start = 1,end = 3))%>%
  group_by(DTA_C)%>% summarise() %>%
  ungroup() %>% mutate(CCanton=as.numeric(DTA_C))

cantones_st
head(cantones_st)

canton_estacional <- data %>% dplyr::select(Canton,CCanton) %>% unique()

G1 <- c('Alajuela','Alajuelita', 'Desamparados', 'San Jose', 'Santa Ana','Turrialba')
G2 <- c('Atenas','Garabito','Orotina')
G3 <- c('Cañas', 'Carillo', 'La Cruz', 'Liberia')
G4 <- c('Esparza','Puntarenas')
G5 <- c('Guacimo','Limon','Matina', 'Pococí', 'Sarapiquí', 'Siquirres')
G6 <- c('SantaCruz','Nicoya')
G7 <- c('Corrredores','Osa')

canton_estacional <- canton_estacional %>% 
  mutate(estacion = case_when( Canton %in% G1 ~ '1',
                               Canton %in% G2 ~ '2',
                               Canton %in% G3 ~ '3',
                               Canton %in% G4 ~ '4',
                               Canton %in% G5 ~ '5',
                               Canton %in% G6 ~ '6',
                               Canton %in% G7 ~ '7',
                               TRUE ~ '8'))


space <- left_join(cantones_st, canton_estacional, by = c("CCanton" = "CCanton"))


rosa          = "#db7590"
cafe          = "#e5b25f"
verdecafe     = "#c6c679"
verdelimon    = "#7ca43f"
verde         = "#2c8e22"
azulTurquesa  = "#3fabc9"
azul          = "#3584da"
morado        = "#998ce1"
purpura       = "#cc4ec0"
gris          = "#636363"

mycolor <- c(rosa,cafe,verdecafe,purpura,verde,morado,azulTurquesa,gris)



month.map <- ggplot(data   = space) +
  geom_sf(mapping     = aes(fill=as.factor(estacion),geometry=geometry),colour="black",
          show.legend = FALSE,size=0.1,alpha=0.8)+
  theme_void()+
  theme(axis.ticks  = element_blank(),axis.text = element_blank())+
  scale_fill_manual(values=mycolor,na.value = 'transparent')


ggsave(month.map, 
       filename = paste0("results_out3mes_final/month.map.jpg"), 
       height = 20, width = 20,bg="white")





