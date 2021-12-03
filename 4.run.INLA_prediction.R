library(tidyverse)
library(timetk)
library(lubridate)
library(gamlss)
library(zoo)
library(boot)
library(dlnm)
library(corrplot)
library(INLA)
library(data.table)


# 0. Dataset preparation -------------------------------------------------

load('./Data/datos_totales.RData')

#introduce clusters
load("./Data/clustering.Rdata")

# Select which cluster to use (clustering by logRR or clustering by logPrecip)
# results_clusters<-results_logRR %>% mutate(cluster=cluster6_logRR) %>%
#                       dplyr::select(Canton,cluster)

results_clusters<-data.frame(Canton= results_logRR$Canton,
                             clusterRR = results_logRR$cluster6_logRR,
                             clusterPrecip = results_logPrecip$cluster6_logPrecip,
                             clusterMVC = results_mvc$cluster6_mvc)

write_csv(results_clusters,file = 'results_clusters.csv')

datos_totales <- datos_totales %>% left_join(results_clusters,by="Canton")

names(datos_totales)
colSums(is.na(datos_totales))
table(datos_totales$Year,datos_totales$Month)

data<- datos_totales %>% filter(Year>=2001 & Year<=2019)
colSums(is.na(data))
# set maximum and minimum lag
maxlag = 12
minlag = 3

#RR
lag_RR <- tsModel::Lag(data$RR, group = data$Canton, k = minlag:maxlag)

# Precipitation
lag_Precip <- tsModel::Lag(data$Precip_t, group = data$Canton, k = minlag:maxlag)
# Nino3SSTA
lag_Nino3SSTA <- tsModel::Lag(data$Nino3SSTA, group = data$Canton, k = minlag:maxlag)
# TNA
lag_TNA <- tsModel::Lag(data$TNA, group = data$Canton, k = minlag:maxlag)

lag_RR<-lag_RR[data$Year>2001,]
lag_Precip<-lag_Precip[data$Year>2001,]
lag_Nino3SSTA<-lag_Nino3SSTA[data$Year>2001,]
lag_TNA<-lag_TNA[data$Year>2001,]

#define data from 2002-2009
data<-data[data$Year>2001,]

# Define time indicator to set 1 to Jan 2002
data <- data %>% group_by(Canton) %>% mutate(time = 1:n()) %>% ungroup()
# total number of months
ntime <- length(unique(data$Month))
# total number of years
nyear <- length(unique(data$Year))
# total number of microregions
nmicro <- length(unique(data$Canton))
# total number of climate clusters
nstate <- length(unique(data$clusterMVC))


# define cross-basis matrix (combining nonlinear exposure and lag functions)
# set lag knots
lagknot = equalknots(minlag:maxlag, 2)

var <- lag_RR
basis_RR <- crossbasis(var,
                       argvar = list(fun = "ns", knots = equalknots(datos_totales$RR, 2)),
                       arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_Precip
basis_Precip <- crossbasis(var,
                           argvar = list(fun = "ns", knots = equalknots(datos_totales$Precip_t, 2)),
                           arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_Nino3SSTA
basis_Nino3SSTA <- crossbasis(var,
                              argvar = list(fun = "ns", knots = equalknots(datos_totales$Nino3SSTA, 2)),
                              arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_TNA
basis_TNA <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(datos_totales$TNA, 2)),
                        arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

colnames(basis_RR) = paste0("b_rr", colnames(basis_RR))
colnames(basis_Precip) = paste0("b_prec", colnames(basis_Precip))
colnames(basis_Nino3SSTA) = paste0("b_nino", colnames(basis_Nino3SSTA))
colnames(basis_TNA) = paste0("b_TNA", colnames(basis_TNA))

data$Year<-data$Year-2001
data$nCanton=as.numeric(factor(data$Canton))

df <- data %>% 
  dplyr::select(Canton,nCanton,time,Year,Month,RR,Cases,Poblacion,OFF,clusterRR,clusterPrecip,clusterMVC) %>%
  mutate(Y = Cases) #%>% 
# bind_cols(df.basis_Precip,df.basis_Nino3SSTA,df.basis_TNA)

df$clusterRR<-as.numeric(df$clusterRR)
df$clusterPrecip<-as.numeric(df$clusterPrecip)
df$clusterMVC<-as.numeric(df$clusterMVC)

df1<-df
df1$Cases[df1$Year==18]<-NA

# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# output plot 
indices_cantones <- c(5,6,7,14,19,23,25,26,27,31)


# 1. Estimation --------------------------------------------------------------

# Models
baseformula1 <- Y ~ 1 + f(Month, replicate = clusterRR, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

baseformula2 <- Y ~ 1 + f(Month, replicate = clusterPrecip, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

baseformula3 <- Y ~ 1 + f(Month, replicate = clusterMVC, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

baseformula4 <- Y ~ 1 + f(Month, replicate = nCanton, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA


formula<-baseformula1
model1 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1 <- inla.rerun(model1)

formula<-baseformula2
model2 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
               control.inla = list(strategy = 'adaptive'), 
               control.compute = list(dic = TRUE, config = FALSE, 
                                      cpo = TRUE, return.marginals = FALSE),
               control.fixed = list(correlation.matrix = TRUE, 
                                    prec.intercept = 1, prec = 1),
               control.predictor = list(link = 1, compute = TRUE), 
               verbose = FALSE)
model2 <- inla.rerun(model2)

formula<-baseformula3
model3 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
               control.inla = list(strategy = 'adaptive'), 
               control.compute = list(dic = TRUE, config = FALSE, 
                                      cpo = TRUE, return.marginals = FALSE),
               control.fixed = list(correlation.matrix = TRUE, 
                                    prec.intercept = 1, prec = 1),
               control.predictor = list(link = 1, compute = TRUE), 
               verbose = FALSE)
model3 <- inla.rerun(model3)

formula<-baseformula4
model4 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
               control.inla = list(strategy = 'adaptive'), 
               control.compute = list(dic = TRUE, config = FALSE, 
                                      cpo = TRUE, return.marginals = FALSE),
               control.fixed = list(correlation.matrix = TRUE, 
                                    prec.intercept = 1, prec = 1),
               control.predictor = list(link = 1, compute = TRUE), 
               verbose = FALSE)
model4 <- inla.rerun(model4)

#model 1: Month |clusterRR
#model 2: Month |clusterPrecip
#model 3: Month |clusterMVC
#model 4: Month |Canton

#CPO
sum(-log(model1$cpo$cpo))
sum(-log(model2$cpo$cpo))
sum(-log(model3$cpo$cpo))
sum(-log(model4$cpo$cpo))

#DIC
model1$dic$dic
model2$dic$dic
model3$dic$dic
model4$dic$dic



# output Model 3 ----------------------------------------------------------

model<-model3
summary(model)


month_effects <- data.table(cbind(rep(unique(df$clusterMVC), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("state_code", "Month")

# plot monthly random effects per state
month.effect.cluster<-month_effects %>% 
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

ggsave(month.effect.cluster, filename = "figs/Model3.month.effect.cluster.jpg", height = 20, width = 20)


canton_effects <- data.frame(canton=unique(df$Canton),
                             effect=model$summary.random$nCanton)

head(canton_effects)

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

ggsave(canton.effect, filename = "figs/Model3.canton.effect.jpg", height = 20, width = 20)


fitted.Case<-model$summary.fitted.values %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

df.final.in <- df.final %>% filter(Year!=18)
df.final.out <- df.final %>% filter(Year==18)


prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = Cases), col = "red") +
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.in, filename = "figs/Model3.prediction.Case.in.jpg", height = 20, width = 20)


prediction.out<-df.final.out %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = Cases), col = "red") +
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.out, filename = "figs/Model3.prediction.Case.out.jpg", height = 20, width = 20)



fitted.RR<-model$summary.linear.predictor %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)


#Recover RR estimates from offset
df.final<-df %>% bind_cols(fitted.RR) %>% 
  mutate(mean=exp(mean-OFF),
         `0.025quant`=exp(`0.025quant`-OFF),
         `0.975quant`=exp(`0.975quant`-OFF)
  )

df.final.in <- df.final %>% filter(Year!=18)
df.final.out <- df.final %>% filter(Year==18)

prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  facet_wrap( ~Canton, scales = "free")

ggsave(prediction.in, filename = "figs/Model3.prediction.RR.in.jpg", height = 20, width = 20)


###

prediction.out<-df.final.out %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.out, filename = "figs/Model3.prediction.RR.out.jpg", height = 20, width = 20)



base.in<- df.final.in %>% 
  mutate(fechas = make_date(Year, Month, 1))%>% 
  dplyr::select(mean,"0.025quant", "0.975quant",RR, fechas, Canton) %>%
  rename(fit = mean,low="0.025quant",up="0.975quant")
base.out<- df.final.out %>% 
  mutate(fechas = make_date(Year, Month, 1))%>% 
  dplyr::select(mean,"0.025quant", "0.975quant",RR, fechas, Canton) %>%
  rename(fit = mean,low="0.025quant",up="0.975quant")

metricas <- function(tabla){
  MSE <- sum((tabla$fit-tabla$RR)^2)
  IS_95 <- mean((tabla$up-tabla$low)+
                  (2/0.05)*(tabla$low-tabla$RR)*(tabla$RR<tabla$low)+
                  (2/0.05)*(tabla$RR-tabla$up)*(tabla$RR>tabla$up))
  return(data.frame(MSE,IS_95))
}

base_grafico_in_g <- base.in %>% group_by(Canton)
base_grafico_out_g <- base.out %>% group_by(Canton)

metricas_tot_in <- base_grafico_in_g %>% group_modify(~metricas(.x))
metricas_tot_out <- base_grafico_out_g %>% group_modify(~metricas(.x))
write_csv(metricas_tot_in,file = 'figs/metricas_clusterMVC_INLA_in.csv')
write_csv(metricas_tot_out,file = 'figs/metricas_clusterMVC_out.csv')



# output Model 4 ----------------------------------------------------------

model<-model4

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

ggsave(month.effect.canton, filename = "figs/Model4.month.effect.canton.jpg", height = 20, width = 20)


canton_effects <- data.frame(canton=unique(df$Canton),
                             effect=model$summary.random$nCanton)

head(canton_effects)

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

ggsave(canton.effect, filename = "figs/Model4.canton.effect.jpg", height = 20, width = 20)



fitted.Case<-model$summary.fitted.values %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

df.final.in <- df.final %>% filter(Year!=18)
df.final.out <- df.final %>% filter(Year==18)


prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = Cases), col = "red") +
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.in, filename = "figs/Model4.prediction.Case.in.jpg", height = 20, width = 20)


prediction.out<-df.final.out %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = Cases), col = "red") +
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.out, filename = "figs/Model4.prediction.Case.out.jpg", height = 20, width = 20)



fitted.RR<-model$summary.linear.predictor %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

#Recover RR estimates from offset
df.final<-df %>% bind_cols(fitted.RR) %>% 
  mutate(mean=exp(mean-OFF),
         `0.025quant`=exp(`0.025quant`-OFF),
         `0.975quant`=exp(`0.975quant`-OFF)
  )

df.final.in <- df.final %>% filter(Year!=18)
df.final.out <- df.final %>% filter(Year==18)

prediction.in <- df.final.in %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  facet_wrap( ~Canton, scales = "free")

ggsave(prediction.in, filename = "figs/Model4.prediction.RR.in.jpg", height = 20, width = 20)


###

prediction.out<-df.final.out %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

ggsave(prediction.out, filename = "figs/Model4.prediction.RR.out.jpg", height = 20, width = 20)



base.in<- df.final.in %>% 
  mutate(fechas = make_date(Year, Month, 1))%>% 
  dplyr::select(mean,"0.025quant", "0.975quant",RR, fechas, Canton) %>%
  rename(fit = mean,low="0.025quant",up="0.975quant")
base.out<- df.final.out %>% 
  mutate(fechas = make_date(Year, Month, 1))%>% 
  dplyr::select(mean,"0.025quant", "0.975quant",RR, fechas, Canton) %>%
  rename(fit = mean,low="0.025quant",up="0.975quant")

metricas <- function(tabla){
  MSE <- sum((tabla$fit-tabla$RR)^2)
  IS_95 <- mean((tabla$up-tabla$low)+
                  (2/0.05)*(tabla$low-tabla$RR)*(tabla$RR<tabla$low)+
                  (2/0.05)*(tabla$RR-tabla$up)*(tabla$RR>tabla$up))
  return(data.frame(MSE,IS_95))
}

base_grafico_in_g <- base.in %>% group_by(Canton)
base_grafico_out_g <- base.out %>% group_by(Canton)

metricas_tot_in <- base_grafico_in_g %>% group_modify(~metricas(.x))
metricas_tot_out <- base_grafico_out_g %>% group_modify(~metricas(.x))
write_csv(metricas_tot_in,file = 'figs/metricas_Canton_INLA_in.csv')
write_csv(metricas_tot_out,file = 'figs/metricas_Canton_INLA_out.csv')



# cases -------------------------------------------------------------------

