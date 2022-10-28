#output

lab <- c("formula.iid", "formula.1.1", "formula.1.2", "formula.1.3", "formula.1.4",
         "formula.2.1", "formula.2.2", "formula.2.3", "formula.2.4",
         "formula.3.1", "formula.3.2", "formula.3.3", "formula.3.4", 
         "formula.1.1a", "formula.1.2a", "formula.1.3a", "formula.1.4a",
         "formula.2.1a", "formula.2.2a", "formula.2.3a", "formula.2.4a",
         "formula.3.1a", "formula.3.2a", "formula.3.3a", "formula.3.4a")


table1 <- data.table(Model = lab, 
                     DIC = NA,
                     logscore = NA)

# create loop to read in model object and extract DIC and CV log score
for (i in 1:length(lab)){
  load(paste0("output_nolineal/", lab[i],".RData"))
  # add model fit statistics from each model to table
  table1$DIC[i] <- round(model$dic$dic,2)
  table1$logscore[i] <- round(-mean(log(model$cpo$cpo), na.rm = T), 4)
}
table1

table1 %>% arrange(DIC)


fwrite(table1, file = "results_nolineal/table_fitting.csv", quote = FALSE, 
       row.names = FALSE) 

# output ------------------------------------------------------------------


# Model 11.1 output ----------------------------------------------------------

source("4.1.load.packages_nolineal.data.R")


selected.model <- "formula.iid"
selected.model <- "formula.1.3a"
selected.model <- "formula.2.3a"
selected.model <- "formula.3.3a"

load(paste0("output_nolineal/",selected.model,".RData"))



month_effects <- data.table(cbind(rep(unique(df$Canton), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("state_code", "Month")

# plot monthly random effects per state
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
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~state_code)

ggsave(month.effect.canton, 
       filename = paste0("results_nolineal/",selected.model,".month.effect.canton.jpg"), 
       height = 20, width = 20)


if(selected.model== "formula.iid"){
  canton_effects <- data.frame(canton=unique(df$Canton),
                               effect=model$summary.random$nCanton)
  
  # plot canton random effects
  
  canton.effect<-canton_effects %>% ggplot(aes(y=canton,x=effect.mean))+
    geom_point()+
    geom_pointrange(aes(xmin = `effect.0.025quant`, xmax = `effect.0.975quant`))+
    theme_bw()+
    labs(x="Efecto aleatorio",x="Canton") +
    geom_vline(xintercept = 0, linetype="dashed", colour="blue")+
    geom_vline(xintercept = -3, linetype="dashed", colour="red")+
    geom_vline(xintercept = 3, linetype="dashed", colour="red")+
    scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,1))
  
}else{
  
  space <- data.table(model$summary.random$nCanton)
  space$year <- rep(min(data$Year):max(data$Year), each = 2*nmicro)
  space$re <- rep(c(rep(1,nmicro),rep(2,nmicro)),nyear)
  space <- space[space$re == 1,]
  space$micro_code <- rep(unique(data$CCanton), nyear)
  mn <-min(space$mean)
  mx <-max(space$mean)
  
  
  ##################
  
  # Add the map geometry to the space dataframe
  space <- left_join(cantones_st, space, by = c("CCanton" = "micro_code"))
  canton.effect <- ggplot() + 
    geom_sf(data = space, aes(fill = mean,geometry=geometry), color = "gray") +
    scale_fill_distiller(palette = "Spectral", direction = -1, 
                         limits = c(min(mn,-mx),max(mx,-mn))) +
    labs(fill = "Contribution to \n log(DIR)") +
    theme_void() +
    facet_wrap(~space$year, ncol = 5)
  
}

#####

ggsave(canton.effect, 
       filename = paste0("results_nolineal/",selected.model,".canton.effect.jpg"), 
       height = 20, width = 20,bg="white")







fitted.Case<-model$summary.fitted.values %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

df.final.in <- df.final %>% filter(nYear!=18)
df.final.out <- df.final %>% filter(nYear==18)


prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) %>%
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
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.in, 
       filename = paste0("results_nolineal/",selected.model,".prediction.Case.in.jpg"), 
       height = 20, width = 20)


prediction.out<-df.final.out %>% filter(nCanton %in% indices_cantones) %>%
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
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.out, 
       filename = paste0("results_nolineal/",selected.model,".prediction.Case.out.jpg"), 
       height = 20, width = 20)



fitted.RR<-model$summary.linear.predictor %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

#Recover RR estimates from offset
df.final<-df %>% bind_cols(fitted.RR) %>% 
  mutate(mean=exp(mean-OFF),
         `0.025quant`=exp(`0.025quant`-OFF),
         `0.975quant`=exp(`0.975quant`-OFF)
  )

df.final.in <- df.final %>% filter(nYear!=18)
df.final.out <- df.final %>% filter(nYear==18)

prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = date, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = date, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = date, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  facet_wrap( ~Canton, scales = "free")

ggsave(prediction.in, 
       filename = paste0("results_nolineal/",selected.model,".prediction.RR.in.jpg"), 
       height = 20, width = 20)


###

prediction.out<-df.final.out %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = date, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = date, y = RR), col = "red") +
  geom_ribbon(aes(x = date, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.out, 
       filename = paste0("results_nolineal/",selected.model,".prediction.RR.out.jpg"), 
       height = 20, width = 20)



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
write_csv(metricas_tot_in,file = paste0("results_nolineal/",selected.model,".metricas.in.csv"))
write_csv(metricas_tot_out,file = paste0("results_nolineal/",selected.model,".metricas.out.csv"))


