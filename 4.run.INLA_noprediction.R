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


# Dataset preparation -------------------------------------------------


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
  dplyr::select(Canton,nCanton,time,Year,Month,RR,Cases,Poblacion,OFF,clusterRR) %>%
  mutate(Y = Cases, cluster=as.numeric(clusterRR)) #%>% 
# bind_cols(df.basis_Precip,df.basis_Nino3SSTA,df.basis_TNA)

# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# output plot 
indices_cantones <- c(5,6,7,14,19,23,25,26,27,31)

# 0. Only temporal effect ----------------------------------------------------

# baseline model
baseformula <- Y ~ 1 + f(Month, replicate = cluster, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) 
# +
#   f(Month, model = "bym2", replicate = Year, graph = "output/map.graph", 
#     scale.model = TRUE, hyper = precision.prior) 

# 
formula0.1 <- update.formula(baseformula, ~. + basis_RR)
formula0.2 <- update.formula(baseformula, ~. + basis_Precip)
formula0.3 <- update.formula(baseformula, ~. + basis_Nino3SSTA)
formula0.4 <- update.formula(baseformula, ~. + basis_TNA)
formula0.5 <- update.formula(baseformula, ~. + basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA)


formula<-baseformula

model0 <- inla(formula = formula, data = df, family = "poisson", offset = OFF,
               control.inla = list(strategy = 'adaptive'), 
               control.compute = list(dic = TRUE, config = FALSE, 
                                      cpo = TRUE, return.marginals = FALSE),
               control.fixed = list(correlation.matrix = TRUE, 
                                    prec.intercept = 1, prec = 1),
               control.predictor = list(link = 1, compute = TRUE), 
               verbose = FALSE)
model0 <- inla.rerun(model0)


formula<-baseformula

model0.0 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model0.0 <- inla.rerun(model0.0)


summary(model0.0)

formula<-formula0.1
model0.1 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model0.1 <- inla.rerun(model0.1)

formula<-formula0.2
model0.2 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model0.2 <- inla.rerun(model0.2)

formula<-formula0.3
model0.3 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model0.3 <- inla.rerun(model0.3)

formula<-formula0.4
model0.4 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model0.4 <- inla.rerun(model0.4)

formula<-formula0.5
model0.5 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model0.5 <- inla.rerun(model0.5)

#CPO
sum(-log(model0$cpo$cpo))
sum(-log(model0.0$cpo$cpo))
sum(-log(model0.1$cpo$cpo))
sum(-log(model0.2$cpo$cpo))
sum(-log(model0.3$cpo$cpo))
sum(-log(model0.4$cpo$cpo))
sum(-log(model0.5$cpo$cpo))

#DIC
model0$dic$dic
model0.0$dic$dic
model0.1$dic$dic
model0.2$dic$dic
model0.3$dic$dic
model0.4$dic$dic
model0.5$dic$dic

plot(model0.5$cpo$cpo)

summary(model0.5)

model0.5$marginals.fixed$b_precv1.l1

model0.5$summary.random$Month

model<-model0.5
month_effects <- data.table(cbind(rep(unique(df$cluster), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("state_code", "Month")

# plot monthly random effects per state
month_effects %>% 
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

model$summary.fixed
model$summary.random
model$summary.linear.predictor

dim(model$summary.fitted.values)
dim(df)

fitted.Case<-model$summary.fitted.values %>% dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

names(df.final)

df.final %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = Cases), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton,scales="free")

fitted.RR<-model$summary.linear.predictor %>% dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.RR)

names(df.final)

df.final %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton)



# 1. Espacial random effect iid cluster----------------------------------------

# baseline model
baseformula <- Y ~ 1 + f(Month, replicate = cluster, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) +
  f(cluster, model="iid")
# +
#   f(Month, model = "bym2", replicate = Year, graph = "output/map.graph", 
#     scale.model = TRUE, hyper = precision.prior) 

# 
formula1.1 <- update.formula(baseformula, ~. + basis_RR)
formula1.2 <- update.formula(baseformula, ~. + basis_Precip)
formula1.3 <- update.formula(baseformula, ~. + basis_Nino3SSTA)
formula1.4 <- update.formula(baseformula, ~. + basis_TNA)
formula1.5 <- update.formula(baseformula, ~. + basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA)


formula<-baseformula

model1.0 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1.0 <- inla.rerun(model1.0)


sum(-log(model0.0$cpo$cpo))
sum(-log(model1.0$cpo$cpo))

model0.0$dic$dic
model1.0$dic$dic

formula<-formula1.1
model1.1 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1.1 <- inla.rerun(model1.1)

formula<-formula1.2
model1.2 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1.2 <- inla.rerun(model1.2)

formula<-formula1.3
model1.3 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1.3 <- inla.rerun(model1.3)

formula<-formula1.4
model1.4 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1.4 <- inla.rerun(model1.4)

formula<-formula1.5
model1.5 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model1.5 <- inla.rerun(model1.5)

#CPO
sum(-log(model1.0$cpo$cpo))
sum(-log(model1.1$cpo$cpo))
sum(-log(model1.2$cpo$cpo))
sum(-log(model1.3$cpo$cpo))
sum(-log(model1.4$cpo$cpo))
sum(-log(model1.5$cpo$cpo))
sum(-log(model0.5$cpo$cpo))

#DIC
model1.0$dic$dic
model1.1$dic$dic
model1.2$dic$dic
model1.3$dic$dic
model1.4$dic$dic
model1.5$dic$dic
model0.5$dic$dic

plot(model1.5$cpo$cpo)

summary(model1.5)

model1.5$marginals.fixed$b_precv1.l1

model1.5$summary.random$Month

model<-model1.5
month_effects <- data.table(cbind(rep(unique(df$cluster), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("state_code", "Month")

# plot monthly random effects per state
month_effects %>% 
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

model$summary.fixed
model$summary.random



cluster_effects <- data.frame(model$summary.random$cluster)

# plot cluster random effects

cluster_effects %>% ggplot(aes(y=ID,x=mean))+
  geom_point()+
  geom_pointrange(aes(xmin = `X0.025quant`, xmax = `X0.975quant`))+
  theme_bw()+
  labs(x="Efecto aleatorio",x="Canton") +
  geom_vline(xintercept = 0, linetype="dashed", colour="blue")+
  geom_vline(xintercept = -3, linetype="dashed", colour="red")+
  geom_vline(xintercept = 3, linetype="dashed", colour="red")+
  scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,1))

df%>% filter(cluster==1)%>% dplyr::select(Canton) %>% unique()


dim(model$summary.fitted.values)


fitted.Case<-model$summary.fitted.values %>% dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

names(df.final)

df.final %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = Cases), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton, scales = "free")


fitted.RR<-model$summary.linear.predictor %>% dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.RR)

#revisar!!!!
df.final1 <- df.final %>% mutate(mean=exp(mean-OFF),
                                 `0.025quant`=exp(`0.025quant`-OFF),
                                 `0.975quant`=exp(`0.975quant`-OFF)
)

df.final<-df.final1


names(df.final)

df.final %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("Cases") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton, scales = "free")



# 2. Espacial random effect iid canton ---------------------------------------

# baseline model
baseformula <- Y ~ 1 + f(Month, replicate = cluster, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) +
  f(nCanton, model="iid")
# +
#   f(Month, model = "bym2", replicate = Year, graph = "output/map.graph", 
#     scale.model = TRUE, hyper = precision.prior) 

# 
formula2.1 <- update.formula(baseformula, ~. + basis_RR)
formula2.2 <- update.formula(baseformula, ~. + basis_Precip)
formula2.3 <- update.formula(baseformula, ~. + basis_Nino3SSTA)
formula2.4 <- update.formula(baseformula, ~. + basis_TNA)
formula2.5 <- update.formula(baseformula, ~. + basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA)



formula<-baseformula

model2.0 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model2.0 <- inla.rerun(model2.0)


formula<-formula2.1
model2.1 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model2.1 <- inla.rerun(model2.1)

formula<-formula2.2
model2.2 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model2.2 <- inla.rerun(model2.2)

formula<-formula2.3
model2.3 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model2.3 <- inla.rerun(model2.3)

formula<-formula2.4
model2.4 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model2.4 <- inla.rerun(model2.4)

formula<-formula2.5
model2.5 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model2.5 <- inla.rerun(model2.5)

#CPO
sum(-log(model2.0$cpo$cpo))
sum(-log(model2.1$cpo$cpo))
sum(-log(model2.2$cpo$cpo))
sum(-log(model2.3$cpo$cpo))
sum(-log(model2.4$cpo$cpo))
sum(-log(model2.5$cpo$cpo))

sum(-log(model0.5$cpo$cpo))
sum(-log(model1.5$cpo$cpo))

#DIC
model2.0$dic$dic
model2.1$dic$dic
model2.2$dic$dic
model2.3$dic$dic
model2.4$dic$dic
model2.5$dic$dic

model0.5$dic$dic
model1.5$dic$dic

plot(model2.5$cpo$cpo)

summary(model2.5)

model2.5$marginals.fixed$b_precv1.l1

model2.5$summary.random$Month

model<-model2.5

model$summary.fixed
model$summary.random
model$summary.linear.predictor

month_effects <- data.table(cbind(rep(unique(df$cluster), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("state_code", "Month")

# plot monthly random effects per state
month_effects %>% 
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



canton_effects <- data.frame(canton=unique(df$Canton),
                                  effect=model$summary.random$nCanton)

head(canton_effects)

# plot canton random effects

  canton_effects %>% ggplot(aes(y=canton,x=effect.mean))+
  geom_point()+
    geom_pointrange(aes(xmin = `effect.0.025quant`, xmax = `effect.0.975quant`))+
  theme_bw()+
  labs(x="Efecto aleatorio",x="Canton") +
  geom_vline(xintercept = 0, linetype="dashed", colour="blue")+
  geom_vline(xintercept = -3, linetype="dashed", colour="red")+
  geom_vline(xintercept = 3, linetype="dashed", colour="red")+
  scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,1))

dim(model$summary.fitted.values)
names(model$summary.fitted.values)

fitted.Case<-model$summary.fitted.values %>% dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

names(df.final)


df.final %>% filter(nCanton %in% indices_cantones) %>%
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

fitted.RR<-model$summary.linear.predictor %>% 
              dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.RR)

df.final1 <- df.final %>% mutate(mean=exp(mean-OFF),
                                 `0.025quant`=exp(`0.025quant`-OFF),
                                 `0.975quant`=exp(`0.975quant`-OFF)
                                 )

df.final<-df.final1

names(df.final)

df.final %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton, scales = "free")




# 3. Espacial random effect iid canton (without temporal effect) ------------------

# baseline model
baseformula <- Y ~ 1 + f(nCanton, model="iid")
# +
#   f(Month, model = "bym2", replicate = Year, graph = "output/map.graph", 
#     scale.model = TRUE, hyper = precision.prior) 

# 
formula3.1 <- update.formula(baseformula, ~. + basis_RR)
formula3.2 <- update.formula(baseformula, ~. + basis_Precip)
formula3.3 <- update.formula(baseformula, ~. + basis_Nino3SSTA)
formula3.4 <- update.formula(baseformula, ~. + basis_TNA)
formula3.5 <- update.formula(baseformula, ~. + basis_RR + basis_Precip + basis_Nino3SSTA + basis_TNA)



formula<-baseformula

model3.0 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model3.0 <- inla.rerun(model3.0)


formula<-formula3.1
model3.1 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model3.1 <- inla.rerun(model3.1)

formula<-formula3.2
model3.2 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model3.2 <- inla.rerun(model3.2)

formula<-formula3.3
model3.3 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model3.3 <- inla.rerun(model3.3)

formula<-formula3.4
model3.4 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model3.4 <- inla.rerun(model3.4)

formula<-formula3.5
model3.5 <- inla(formula = formula, data = df, family = "nbinomial", offset = OFF,
                 control.inla = list(strategy = 'adaptive'), 
                 control.compute = list(dic = TRUE, config = FALSE, 
                                        cpo = TRUE, return.marginals = FALSE),
                 control.fixed = list(correlation.matrix = TRUE, 
                                      prec.intercept = 1, prec = 1),
                 control.predictor = list(link = 1, compute = TRUE), 
                 verbose = FALSE)
model3.5 <- inla.rerun(model3.5)

#CPO
sum(-log(model3.0$cpo$cpo))
sum(-log(model3.1$cpo$cpo))
sum(-log(model3.2$cpo$cpo))
sum(-log(model3.3$cpo$cpo))
sum(-log(model3.4$cpo$cpo))
sum(-log(model3.5$cpo$cpo))

sum(-log(model2.5$cpo$cpo))

#DIC
model3.0$dic$dic
model3.1$dic$dic
model3.2$dic$dic
model3.3$dic$dic
model3.4$dic$dic
model3.5$dic$dic

model0.5$dic$dic
model1.5$dic$dic

plot(model3.5$cpo$cpo)

summary(model3.5)

model3.5$marginals.fixed$b_precv1.l1

model3.5$summary.random$Month

model<-model3.5

model$summary.fixed
model$summary.random
model$summary.linear.predictor

month_effects <- data.table(cbind(rep(unique(df$cluster), each = 12),
                                  model$summary.random$Month))
names(month_effects)[1:2] <- c("state_code", "Month")

# plot monthly random effects per state
month_effects %>% 
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



canton_effects <- data.frame(canton=unique(df$Canton),
                             effect=model$summary.random$nCanton)

head(canton_effects)

# plot canton random effects

canton_effects %>% ggplot(aes(y=canton,x=effect.mean))+
  geom_point()+
  geom_pointrange(aes(xmin = `effect.0.025quant`, xmax = `effect.0.975quant`))+
  theme_bw()+
  labs(x="Efecto aleatorio",x="Canton") +
  geom_vline(xintercept = 0, linetype="dashed", colour="blue")+
  geom_vline(xintercept = -3, linetype="dashed", colour="red")+
  geom_vline(xintercept = 3, linetype="dashed", colour="red")+
  scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,1))

dim(model$summary.fitted.values)
names(model$summary.fitted.values)

fitted.Case<-model$summary.fitted.values %>% dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.Case)

names(df.final)


df.final %>% filter(nCanton %in% indices_cantones) %>%
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

fitted.RR<-model$summary.linear.predictor %>% 
  dplyr::select(mean,`0.025quant`,`0.975quant`)

df.final<-df %>% bind_cols(fitted.RR)

df.final1 <- df.final %>% mutate(mean=exp(mean-OFF),
                                 `0.025quant`=exp(`0.025quant`-OFF),
                                 `0.975quant`=exp(`0.975quant`-OFF)
)

df.final<-df.final1

names(df.final)

df.final %>% filter(nCanton %in% indices_cantones) %>%
  ggplot() + 
  geom_ribbon(aes(x = time, ymin = `0.025quant`, ymax = `0.975quant`), 
              fill = "cadetblue4", alpha = 0.5) + 
  geom_line(aes(x = time, y = `mean`), col = "cadetblue4") +
  geom_line(aes(x = time, y = RR), col = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("time") +
  ylab("RR") +
  theme_bw() + 
  # organise by state name in grid file
  facet_wrap( ~Canton, scales = "free")



