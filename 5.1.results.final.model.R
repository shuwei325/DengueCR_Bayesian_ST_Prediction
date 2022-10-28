source("4.1.load.packages.data_out3mes.R")

selected.model <- "formula.1.2a"  #modicar el mapa de cantones por año

load(paste0("output_out3mes/",selected.model,".RData"))


# Output ------------------------------------------------------------------

month_effects <- data.table(cbind(rep(unique(df$Canton), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("Canton", "Month")


G1 <- c('Alajuela','Alajuelita', 'Desamparados', 'San Jose', 'Santa Ana','Turrialba')
G2 <- c('Atenas','Garabito','Orotina')
G3 <- c('Cañas', 'Carillo', 'La Cruz', 'Liberia')
G4 <- c('Esparza','Puntarenas')
G5 <- c('Guacimo','Limon','Matina', 'Pococí', 'Sarapiquí', 'Siquirres')
G6 <- c('SantaCruz','Nicoya')
G7 <- c('Corrredores','Osa')

month_effects <- month_effects %>% 
  mutate(estacion = case_when( Canton %in% G1 ~ '1',
                               Canton %in% G2 ~ '2',
                               Canton %in% G3 ~ '3',
                               Canton %in% G4 ~ '4',
                               Canton %in% G5 ~ '5',
                               Canton %in% G6 ~ '6',
                               Canton %in% G7 ~ '7',
                               TRUE ~ '8'))

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

hacky_df <- month_effects %>% select(Canton, estacion) %>% unique()

  month.effect.canton<-month_effects %>%
  ggplot() + 
  geom_ribbon(aes(x = Month, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = Month, y = `mean`), col = "cadetblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("Month") +
  ylab("Contribution to log(DIR)") +
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  facet_wrap( ~Canton)+
  coord_cartesian(clip="off", ylim=c(-1, 1)) +
  geom_rect(
    data=hacky_df,
    aes(xmin=-Inf, xmax=Inf,
        ymin=1, ymax=1.5,     # totally defined by trial-and-error
        fill=estacion),alpha=0.4)+
  scale_fill_manual("Groups",values = c("1" = mycolor[1], 
                                        "2" = mycolor[2], 
                                        "3" = mycolor[3], 
                                        "4" = mycolor[4],
                                        "5" = mycolor[5],
                                        "6" = mycolor[6],
                                        "7" = mycolor[7],
                                        "8" = mycolor[8])) +
  theme_bw() + 
  theme( axis.text = element_text(size = 20),
         strip.text = element_text(size = 30),
         strip.background = element_rect(fill=NA),
         legend.title = element_text( size = 30),
         legend.text = element_text(size = 30))


ggsave(month.effect.canton, 
       filename = paste0("results_out3mes_final/month.effect.canton.jpg"), 
       height = 20, width = 20)



# maps --------------------------------------------------------------------

space <- data.table(model$summary.random$nCanton)
space$year <- rep(min(data$Year):max(data$Year), each = nmicro)
space$micro_code <- rep(unique(data$CCanton), nyear)
mn <-min(space$mean)
mx <-max(space$mean)


# Add the map geometry to the space dataframe

codCantones <- cantones_st%>%
  select(CCanton)%>%distinct()%>% 
  slice(rep(1:n(), 20))%>%mutate(Year=c(rep(2002,81),rep(2003,81),
                                        rep(2004,81),rep(2005,81),rep(2006,81),
                                        rep(2007,81),rep(2008,81),rep(2009,81),
                                        rep(2010,81),rep(2011,81),rep(2012,81),
                                        rep(2013,81),rep(2014,81),rep(2015,81),
                                        rep(2016,81),rep(2017,81),rep(2018,81),
                                        rep(2019,81),rep(2020,81),rep(2021,81)))

space <- full_join(codCantones, space, by = c("CCanton" = "micro_code","Year" = "year")) 

canton.effect <- ggplot() + 
  geom_sf(data = space, aes(fill = mean,geometry=geometry), color = "black") +
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value = 'transparent', 
                       limits = c(min(mn,-mx),max(mx,-mn))) +
  labs(fill = "Spatial \n random effect") +
  theme_void() +
  facet_wrap(~space$Year, ncol = 5)+
  theme( legend.key.size = unit(1, 'cm'),
         legend.title = element_text(size=30),
         legend.text = element_text(size=30),
         strip.text = element_text(size = 30))

ggsave(canton.effect, 
       filename = paste0("results_out3mes_final/canton.effect.jpg"), 
       height = 20, width = 20,bg="white")



# Prediction - Cases --------------------------------------------------------------


fitted.Case<-model$summary.fitted.values %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

df.final.in <- df.final %>% filter(nYear!=20)
df.final.out <- df.final %>% filter(nYear==20)


cantones <- df.final.in %>% select(Canton,nCanton) %>% unique()
# output plot 


indices_cantones <- c(30,21,18,1,14,23,24, 16, 13) #3 best and 3 worst


prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) 

prediction.in$Canton <- factor(prediction.in$Canton,
                           levels = c("Talamanca", "Perez Zeledón", "Orotina",
                                      "Alajuela", "Limon", "Puntarenas",
                                      "Quepos", "Montes de Oro","Liberia"))
prediction.in<-prediction.in %>%
  ggplot() + 
  geom_line(aes(x = date, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = date, y = Cases), col = "red") +
  geom_ribbon(aes(x = date, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free") +
  theme( axis.text = element_text( size = 20 ),
         axis.title = element_text( size = 30, face = "bold" ),
         legend.position="none",
         # The new stuff
         strip.text = element_text(size = 40))


ggsave(prediction.in, 
       filename = paste0("results_out3mes_final/prediction.Case.in.jpg"), 
       height = 20, width = 20)



prediction.out <- df.final.out %>% filter(nCanton %in% indices_cantones)

prediction.out$Canton <- factor(prediction.out$Canton,
                               levels = c("Talamanca", "Perez Zeledón", "Orotina",
                                          "Alajuela", "Limon", "Puntarenas",
                                          "Quepos", "Montes de Oro","Liberia"))


prediction.out<-prediction.out %>%
  ggplot() + 
  geom_line(aes(x = date, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = date, y = Cases), col = "red") +
  geom_ribbon(aes(x = date, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") + 
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")+
  theme( axis.text = element_text( size = 20 ),
         axis.title = element_text( size = 30, face = "bold" ),
         legend.position="none",
         # The new stuff
         strip.text = element_text(size = 40))

ggsave(prediction.out, 
       filename = paste0("results_out3mes_final/prediction.Case.out.jpg"), 
       height = 20, width = 20)





# Prediction - RR --------------------------------------------------------------


fitted.RR<-model$summary.linear.predictor %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

#Recover RR estimates from offset
df.final<-df %>% bind_cols(fitted.RR) %>% 
  mutate(mean=exp(mean-OFF),
         `0.025quant`=exp(`0.025quant`-OFF),
         `0.975quant`=exp(`0.975quant`-OFF)
  )

df.final <- df.final %>% mutate(dif = abs(RR-mean)/RR,
                                ampl.interv = `0.975quant` - `0.025quant`)


df.final.in <- df.final %>% filter(nYear!=20)
df.final.out <- df.final %>% filter(nYear==20)


prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) 

prediction.in$Canton <- factor(prediction.in$Canton,
                               levels = c("Talamanca", "Perez Zeledón", "Orotina",
                                          "Alajuela", "Limon", "Puntarenas",
                                          "Quepos", "Montes de Oro","Liberia"))
prediction.in<-prediction.in %>%
  ggplot() + 
  geom_ribbon(aes(x = date, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = date, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = date, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  facet_wrap( ~Canton, scales = "free")+ 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")+
  theme( axis.text = element_text( size = 20 ),
         axis.title = element_text( size = 30, face = "bold" ),
         legend.position="none",
         # The new stuff
         strip.text = element_text(size = 40))

ggsave(prediction.in, 
       filename = paste0("results_out3mes_final/prediction.RR.in.jpg"), 
       height = 20, width = 20)


###

prediction.out <- df.final.out %>% filter(nCanton %in% indices_cantones)

prediction.out$Canton <- factor(prediction.out$Canton,
                                levels = c("Talamanca", "Perez Zeledón", "Orotina",
                                           "Alajuela", "Limon", "Puntarenas",
                                           "Quepos", "Montes de Oro","Liberia"))


prediction.out<-prediction.out %>%
  ggplot() + 
  geom_line(aes(x = date, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = date, y = RR), col = "red") +
  geom_ribbon(aes(x = date, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")+ 
  theme( axis.text = element_text( size = 20 ),
         axis.title = element_text( size = 30, face = "bold" ),
         legend.position="none",
         # The new stuff
         strip.text = element_text(size = 40))

ggsave(prediction.out, 
       filename = paste0("results_out3mes_final/prediction.RR.out.jpg"), 
       height = 20, width = 20)




# Prediction map RR ----------------------------------------------------------

# in

j <- 2002:2020

for (i in 11:19){
#for (i in seq_along(j)){
year.i <- j[i]   # 2017, 2018, 2019, 2020

df.final.in <- df.final %>% filter(Year%in%c(year.i))
mn <-min(df.final.in$mean)
mx <-max(df.final.in$mean)

codCantones <- cantones_st%>%
  select(CCanton)%>%distinct()%>% 
  slice(rep(1:n(), 12))%>%mutate(Month = rep(1:12, each=81))


RR.in <- full_join(codCantones, df.final.in, by = c("CCanton" = "CCanton","Month" = "Month")) 

prediction.map.in <- ggplot() + 
  geom_sf(data = RR.in, aes(fill = mean,geometry=geometry), color = "black") +
  scale_fill_gradient(low = "cyan", high = "red", na.value = 'transparent'
                      ,limits = c(mn,mx))+
  labs(fill = "RR mean") +
  #ggtitle(as.character(year.i)) +
  theme_void() +
  facet_wrap(~RR.in$Month, ncol = 3) +
  theme( legend.key.size = unit(1, 'cm'),
         legend.title = element_text(size=30),
         legend.text = element_text(size=30),
         strip.text = element_text(size = 30))

ggsave(prediction.map.in, 
       filename = paste0("results_out3mes_final/prediction.map.RR.in.",year.i,".jpg"), 
       height = 20, width = 20)



dif.map.in <- ggplot() + 
  geom_sf(data = RR.in, aes(fill = dif,geometry=geometry), color = "black") +
  scale_fill_gradient(low = "cyan", high = "red", na.value = 'transparent')+
  labs(fill = "Absolute \n percentage error") +
  ggtitle(as.character(year.i)) +
  theme_void() +
  facet_wrap(~RR.in$Month, ncol = 3) +
  theme(  plot.title = element_text(size = 40, face = "bold"),
          legend.key.size = unit(1, 'cm'),
         legend.title = element_text(size=30),
         legend.text = element_text(size=30),
         strip.text = element_text(size = 30))

ggsave(dif.map.in, 
       filename = paste0("results_out3mes_final/dif.map.RR.in.",year.i,".jpg"), 
       height = 20, width = 20)
}



# out


year.i <- 2021


df.final.out <- df.final %>% filter(Year%in%c(year.i))
mn <-min(df.final.out$mean)
mx <-max(df.final.out$mean)


codCantones <- cantones_st%>%
  select(CCanton)%>%distinct()%>% 
  slice(rep(1:n(), 3))%>%mutate(Month=rep(1:3, each=81))

RR.out <- full_join(codCantones, df.final.out, by = c("CCanton" = "CCanton","Month" = "Month")) 

prediction.map.out <- ggplot() + 
  geom_sf(data = RR.out, aes(fill = mean,geometry=geometry), color = "black") +
  scale_fill_gradient(low = "cyan", high = "red", na.value = 'transparent', 
                      limits = c(mn,mx))+
  labs(fill = "RR") +
  theme_void() +
  facet_wrap(~RR.out$Month, ncol = 2) +
  theme( legend.key.size = unit(1, 'cm'),
         legend.title = element_text(size=30),
         legend.text = element_text(size=30),
         strip.text = element_text(size = 30))


ggsave(prediction.map.out, 
       filename = paste0("results_out3mes_final/prediction.map.RR.out.jpg"), 
       height = 20, width = 20)



dif.map.out <- ggplot() + 
  geom_sf(data = RR.out, aes(fill = dif,geometry=geometry), color = "black") +
  scale_fill_gradient(low = "cyan", high = "red", na.value = 'transparent')+
  labs(fill = "Absolute \n percentage error") +
  theme_void() +
  facet_wrap(~RR.out$Month, ncol = 2) +
  theme( legend.key.size = unit(1, 'cm'),
         legend.title = element_text(size=30),
         legend.text = element_text(size=30),
         strip.text = element_text(size = 30))


ggsave(dif.map.out, 
       filename = paste0("results_out3mes_final/dif.map.RR.out.jpg"), 
       height = 20, width = 20)





# Métricas ----------------------------------------------------------------


base.in<- df.final.in %>% 
  mutate(fechas = make_date(Year, Month, 1))%>% 
  dplyr::select(mean,"0.025quant", "0.975quant",RR, fechas, Canton) %>%
  rename(fit = mean,low="0.025quant",up="0.975quant")
base.out<- df.final.out %>% 
  mutate(fechas = make_date(Year, Month, 1))%>% 
  dplyr::select(mean,"0.025quant", "0.975quant",RR, fechas, Canton) %>%
  rename(fit = mean,low="0.025quant",up="0.975quant")

metricas <- function(tabla){
  MSE <- sum((tabla$fit-tabla$RR)^2)/mean(tabla$RR)
  IS_95 <- mean((tabla$up-tabla$low)+
                  (2/0.05)*(tabla$low-tabla$RR)*(tabla$RR<tabla$low)+
                  (2/0.05)*(tabla$RR-tabla$up)*(tabla$RR>tabla$up))/mean(tabla$RR)
  return(data.frame(MSE,IS_95))
}

base_grafico_in_g <- base.in %>% group_by(Canton)
base_grafico_out_g <- base.out %>% group_by(Canton)

metricas_tot_in <- base_grafico_in_g %>% group_modify(~metricas(.x))
metricas_tot_out <- base_grafico_out_g %>% group_modify(~metricas(.x))
write_csv(metricas_tot_in,file = paste0("results_out3mes/",selected.model,".metricas.in.csv"))
write_csv(metricas_tot_out,file = paste0("results_out3mes/",selected.model,".metricas.out.csv"))













