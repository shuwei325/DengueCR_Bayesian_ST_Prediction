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
                             clusterMVC3 = results_mvc$cluster3_mvc,
                             clusterMVC4 = results_mvc$cluster4_mvc,
                             clusterMVC5 = results_mvc$cluster5_mvc,
                             clusterMVC6 = results_mvc$cluster6_mvc)

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

###########################
# Nino12SSTA
lag_Nino12SSTA <- tsModel::Lag(data$Nino12SSTA, group = data$Canton, k = minlag:maxlag)
# Nino4SSTA
lag_Nino4SSTA <- tsModel::Lag(data$Nino4SSTA, group = data$Canton, k = minlag:maxlag)
# Nino34SSTA
lag_Nino34SSTA <- tsModel::Lag(data$Nino34SSTA, group = data$Canton, k = minlag:maxlag)
# EVI
lag_EVI <- tsModel::Lag(data$EVI, group = data$Canton, k = minlag:maxlag)
# NDVI
lag_NDVI <- tsModel::Lag(data$NDVI, group = data$Canton, k = minlag:maxlag)
# NDWI
lag_NDWI <- tsModel::Lag(data$NDWI, group = data$Canton, k = minlag:maxlag)
# LSD
lag_LSD <- tsModel::Lag(data$LSD, group = data$Canton, k = minlag:maxlag)
# LSN
lag_LSN <- tsModel::Lag(data$LSN, group = data$Canton, k = minlag:maxlag)
######################



lag_RR<-lag_RR[data$Year>2001,]
lag_Precip<-lag_Precip[data$Year>2001,]
lag_Nino3SSTA<-lag_Nino3SSTA[data$Year>2001,]
lag_TNA<-lag_TNA[data$Year>2001,]

##########################
lag_Nino12SSTA <- lag_Nino12SSTA[data$Year>2001,]
lag_Nino4SSTA <- lag_Nino4SSTA[data$Year>2001,]
lag_Nino34SSTA <- lag_Nino34SSTA[data$Year>2001,]
lag_EVI <- lag_EVI[data$Year>2001,]
lag_NDVI <- lag_NDVI[data$Year>2001,]
lag_NDWI <- lag_NDWI[data$Year>2001,]
lag_LSD <- lag_LSD[data$Year>2001,]
lag_LSN <- lag_LSN[data$Year>2001,]
##########################



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
colnames(basis_Nino3SSTA) = paste0("b_nino3", colnames(basis_Nino3SSTA))
colnames(basis_TNA) = paste0("b_TNA", colnames(basis_TNA))

##########################
var <- lag_Nino12SSTA
basis_Nino12SSTA <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(datos_totales$Nino12SSTA, 2)),
                        arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_Nino4SSTA
basis_Nino4SSTA <- crossbasis(var,
                               argvar = list(fun = "ns", knots = equalknots(datos_totales$Nino4SSTA, 2)),
                               arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_Nino34SSTA
basis_Nino34SSTA <- crossbasis(var,
                               argvar = list(fun = "ns", knots = equalknots(datos_totales$Nino34SSTA, 2)),
                               arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_EVI
basis_EVI <- crossbasis(var,
                               argvar = list(fun = "ns", knots = equalknots(datos_totales$EVI, 2)),
                               arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_NDVI
basis_NDVI <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(datos_totales$NDVI, 2)),
                        arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_NDWI
basis_NDWI <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(datos_totales$NDWI, 2)),
                        arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_LSD
basis_LSD <- crossbasis(var,
                         argvar = list(fun = "ns", knots = equalknots(datos_totales$LSD, 2)),
                         arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

var <- lag_LSN
basis_LSN <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(datos_totales$LSN, 2)),
                        arglag = list(fun = "ns", knots = (maxlag-minlag)/2))

colnames(basis_Nino12SSTA) = paste0("b_nino12", colnames(basis_Nino12SSTA))
colnames(basis_Nino4SSTA) = paste0("b_nino3", colnames(basis_Nino4SSTA))
colnames(basis_Nino34SSTA) = paste0("b_nino34", colnames(basis_Nino34SSTA))
colnames(basis_EVI) = paste0("b_EVI", colnames(basis_EVI))
colnames(basis_NDVI) = paste0("b_NDVI", colnames(basis_NDVI))
colnames(basis_NDWI) = paste0("b_NDWI", colnames(basis_NDWI))
colnames(basis_LSD) = paste0("b_LSD", colnames(basis_LSD))
colnames(basis_LSN) = paste0("b_LSN", colnames(basis_LSN))

#######################################



data$Year<-data$Year-2001
data$nCanton=as.numeric(factor(data$Canton))

df <- data %>% 
  dplyr::select(Canton,nCanton,time,Year,Month,RR,Cases,Poblacion,OFF,clusterRR,clusterPrecip,clusterMVC3,clusterMVC4,clusterMVC5,clusterMVC6) %>%
  mutate(Y = Cases) 


df$clusterRR<-as.numeric(df$clusterRR)
df$clusterPrecip<-as.numeric(df$clusterPrecip)
df$clusterMVC3<-as.numeric(df$clusterMVC3)
df$clusterMVC4<-as.numeric(df$clusterMVC4)
df$clusterMVC5<-as.numeric(df$clusterMVC5)
df$clusterMVC6<-as.numeric(df$clusterMVC6)

df1<-df
df1$Cases[df1$Year==18]<-NA

# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# output plot 
indices_cantones <- c(5,6,7,14,19,23,25,26,27,31)


# 1. Estimation --------------------------------------------------------------


# include formula and set defaults for data, family (to allow other prob dist models e.g. Poisson) and config (to allow for sampling)
mymodel <- function(formula, data = df1, family = "nbinomial", config = FALSE){
  model <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                control.inla = list(strategy = 'adaptive'), 
                control.compute = list(dic = TRUE, config = FALSE, 
                                       cpo = TRUE, return.marginals = FALSE),
                control.fixed = list(correlation.matrix = TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                verbose = FALSE)
  model <- inla.rerun(model)
  return(model)
}


# Models
# Comparing different cluster and Canton options
formula_clusterRR <- Y ~ 1 + f(Month, replicate = clusterRR, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

formula_clusterPrecip <- Y ~ 1 + f(Month, replicate = clusterPrecip, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

formula_clusterMVC3 <- Y ~ 1 + f(Month, replicate = clusterMVC3, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

formula_clusterMVC4 <- Y ~ 1 + f(Month, replicate = clusterMVC4, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

formula_clusterMVC5 <- Y ~ 1 + f(Month, replicate = clusterMVC5, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

formula_clusterMVC6 <- Y ~ 1 + f(Month, replicate = clusterMVC6, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

formula_canton <- Y ~ 1 + f(Month, replicate = nCanton, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") + 
  basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA

##############

formulas <- list(formula_clusterRR, formula_clusterPrecip, 
      formula_clusterMVC3, formula_clusterMVC4, formula_clusterMVC5, 
      formula_clusterMVC6, formula_canton)
lab <- c("mod_clusterRR", "mod_clusterPrecip", 
         "mod_clusterMVC3", "mod_clusterMVC4", "mod_clusterMVC5", 
         "mod_clusterMVC6", "mod_canton")

K <- length(formulas)
models <- lapply(1:K, 
                 function(i) {
                   model <- mymodel(formulas[[i]], data=df1)
                   save(model, file = paste0("output/", lab[i],".RData"))
                   }
                 )


# create table to store DIC and select best model 
table0 <- data.table(Model  = lab, 
                     DIC = NA, CPO = NA)

for(i in 1:length(formulas)){
  load(paste0("output/",lab[i],".RData"))
  table0$DIC[i] <- round(model$dic$dic)
  table0$CPO[i] <- round(sum(-log(model$cpo$cpo)))
}

table0 %>% arrange(.,DIC)



# 2. Estimation 2 ---------------------------------------------------------


baseformula <- Y ~ 1 + f(Month, replicate = clusterMVC4, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") +
                basis_RR + basis_Precip

formula_Nino3 <- update.formula(baseformula, ~. + basis_Nino3SSTA)
formula_Nino12 <- update.formula(baseformula, ~. + basis_Nino12SSTA)
formula_Nino4 <- update.formula(baseformula, ~. + basis_Nino4SSTA)
formula_Nino34 <- update.formula(baseformula, ~. + basis_Nino34SSTA)

formulas <- list(baseformula, formula_Nino3,
                 formula_Nino12,formula_Nino4,formula_Nino34)
lab <- c("2_base_mod", "2_mod_Nino3",
         "2_mod_Nino12","2_mod_Nino4","2_mod_Nino34")

K <- length(formulas)
models <- lapply(1:K, 
                 function(i) {
                   model <- mymodel(formulas[[i]], data=df1)
                   save(model, file = paste0("output/", lab[i],".RData"))
                 }
)

# create table to store DIC and select best model 
table1 <- data.table(Model  = lab, 
                     DIC = NA, CPO = NA)

for(i in 1:length(formulas)){
  load(paste0("output/",lab[i],".RData"))
  table1$DIC[i] <- round(model$dic$dic)
  table1$CPO[i] <- round(sum(-log(model$cpo$cpo)))
}

table1 %>% arrange(.,DIC)



# 3. Estimation 3 ---------------------------------------------------------


formula_final1 <- Y ~ 1 + f(Month, replicate = clusterMVC4, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") +
  basis_RR + basis_Precip + basis_Nino12SSTA + basis_TNA + 
  basis_NDVI + basis_LSD 

formula_final2 <- Y ~ 1 + f(Month, replicate = nCanton, model = "rw1", cyclic = TRUE, constr = TRUE,
                            scale.model = TRUE,  hyper = precision.prior) + f(nCanton, model="iid") +
  basis_RR + basis_Precip + basis_Nino12SSTA + basis_TNA + 
  basis_NDVI + basis_LSD 

formulas <- list(formula_final1,formula_final2)
lab <- c("final_mod1","final_mod2")


K <- length(formulas)
models <- lapply(1:K, 
                 function(i) {
                   model <- mymodel(formulas[[i]], data=df1)
                   save(model, file = paste0("output/", lab[i],".RData"))
                 }
)


# create table to store DIC and select best model 
table2 <- data.table(Model  = lab, 
                     DIC = NA, CPO = NA)

for(i in 1:length(formulas)){
  load(paste0("output/",lab[i],".RData"))
  table2$DIC[i] <- round(model$dic$dic)
  table2$CPO[i] <- round(sum(-log(model$cpo$cpo)))
}

bind_rows(table0,table1,table2) %>% arrange(.,DIC)

# 4. output Model 1 ----------------------------------------------------------

i <- 1
mod.name <- c("final_mod1","final_mod2")
load(paste0("output/", mod.name[i],".RData"))

summary(model)


month_effects <- data.table(cbind(rep(unique(df1$clusterMVC4), each = 12),
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

ggsave(month.effect.cluster, filename = "results/ModelMVC4.month.effect.cluster.jpg", height = 20, width = 20)


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

ggsave(canton.effect, filename = "results/ModelMVC4.canton.effect.jpg", height = 20, width = 20)


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

ggsave(prediction.in, filename = "results/ModelMVC4.prediction.Case.in.jpg", height = 20, width = 20)


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

ggsave(prediction.out, filename = "results/ModelMVC4.prediction.Case.out.jpg", height = 20, width = 20)



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

ggsave(prediction.in, filename = "results/ModelMVC4.prediction.RR.in.jpg", height = 20, width = 20)


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

ggsave(prediction.out, filename = "results/ModelMVC4.prediction.RR.out.jpg", height = 20, width = 20)



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
write_csv(metricas_tot_in,file = 'results/ModelMVC4_metricas_in.csv')
write_csv(metricas_tot_out,file = 'results/ModelMVC4_metricas_out.csv')

# 4. output Model 2 ----------------------------------------------------------

i <- 2
mod.name <- c("final_mod1","final_mod2")
load(paste0("output/", mod.name[i],".RData"))

summary(model)

month_effects <- data.table(cbind(rep(unique(df1$Canton), each = 12),
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

ggsave(month.effect.canton, filename = "results/ModelCanton.month.effect.canton.jpg", height = 20, width = 20)


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

ggsave(canton.effect, filename = "results/ModelCanton.canton.effect.jpg", height = 20, width = 20)



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

ggsave(prediction.in, filename = "results/ModelCanton.prediction.Case.in.jpg", height = 20, width = 20)


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

ggsave(prediction.out, filename = "results/ModelCanton.prediction.Case.out.jpg", height = 20, width = 20)



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

ggsave(prediction.in, filename = "results/ModelCanton.prediction.RR.in.jpg", height = 20, width = 20)


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

ggsave(prediction.out, filename = "results/ModelCanton.prediction.RR.out.jpg", height = 20, width = 20)



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
write_csv(metricas_tot_in,file = 'results/ModelCanton_Canton_in.csv')
write_csv(metricas_tot_out,file = 'results/ModelCanton_Canton_out.csv')
# cases -------------------------------------------------------------------

