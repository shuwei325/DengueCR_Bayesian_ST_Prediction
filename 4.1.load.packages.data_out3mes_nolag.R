library(tidyverse)
library(timetk)
library(lubridate)
library(zoo)
library(boot)
library(dlnm)
library(corrplot)
library(INLA)
library(data.table)


# 0. Dataset preparation -------------------------------------------------

load('./Data/datos_totales.RData')


#Archivos para el mapa
distritos_st <-st_read('./Data/distritos_shape/Distritos_de_Costa_Rica.shp') %>%
  filter(CODIGO!=60110)

cantones_st <- distritos_st %>%
  mutate(DTA_C = str_sub(CODIGO,start = 1,end = 3))%>%
  group_by(DTA_C)%>% summarise() %>%
  ungroup() %>% mutate(CCanton=as.numeric(DTA_C))



#introduce clusters
load("./Data/clustering.Rdata")

#introduce clusters usando wavelets
load("./Data/wavelets_clusters.Rdata")
colnames(wavelets_clusters)[2:6]<-c("cluster3_wl","cluster4_wl",
                                    "cluster5_wl","cluster6_wl","cluster7_wl")

data_clusters   <- wavelets_clusters %>%
  left_join(results_mvc, by="Canton")

data_clusters <- data_clusters %>% dplyr::select(Canton,cluster6_mvc,cluster6_wl,cluster7_wl)


datos_totales <- datos_totales %>% left_join(data_clusters,by="Canton")


names(datos_totales)
colSums(is.na(datos_totales))
table(datos_totales$Year,datos_totales$Month)

data<- datos_totales %>% filter(Year>=2001 & Year<=2021)
data<- data %>% filter(Year<=2020 | Month<=3 )
colSums(is.na(data))
# set maximum and minimum lag
maxlag = 12
minlag = 0

maxlag.rr = 1
minlag.rr = 1

#RR
lag_RR <- tsModel::Lag(data$RR, group = data$Canton, k = minlag.rr:maxlag.rr)

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



lag_RR<-lag_RR[data$Year>2001]
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
#nstate <- length(unique(data$clusterMVC))


# define cross-basis matrix (combining nonlinear exposure and lag functions)
# set lag knots
lagknot = equalknots(minlag:maxlag, 2)


var <- lag_RR
basis_RR <- crossbasis(var,
                       argvar = list(fun = "lin"),
                       arglag = list(fun = "lin"))

var <- lag_Precip
basis_Precip <- crossbasis(var,
                           argvar = list(fun = "lin"),
                           arglag = list(fun = "lin"))

var <- lag_Nino3SSTA
basis_Nino3SSTA <- crossbasis(var,
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "lin"))

var <- lag_TNA
basis_TNA <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "lin"))

colnames(basis_RR) = paste0("b_rr", colnames(basis_RR))
colnames(basis_Precip) = paste0("b_prec", colnames(basis_Precip))
colnames(basis_Nino3SSTA) = paste0("b_nino3", colnames(basis_Nino3SSTA))
colnames(basis_TNA) = paste0("b_TNA", colnames(basis_TNA))

##########################
var <- lag_Nino12SSTA
basis_Nino12SSTA <- crossbasis(var,
                               argvar = list(fun = "lin"),
                               arglag = list(fun = "lin"))

var <- lag_Nino4SSTA
basis_Nino4SSTA <- crossbasis(var,
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "lin"))

var <- lag_Nino34SSTA
basis_Nino34SSTA <- crossbasis(var,
                               argvar = list(fun = "lin"),
                               arglag = list(fun = "lin"))

var <- lag_EVI
basis_EVI <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "lin"))

var <- lag_NDVI
basis_NDVI <- crossbasis(var,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "lin"))

var <- lag_NDWI
basis_NDWI <- crossbasis(var,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "lin"))

var <- lag_LSD
basis_LSD <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "lin"))

var <- lag_LSN
basis_LSN <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "lin"))

colnames(basis_Nino12SSTA) = paste0("b_nino12", colnames(basis_Nino12SSTA))
colnames(basis_Nino4SSTA) = paste0("b_nino3", colnames(basis_Nino4SSTA))
colnames(basis_Nino34SSTA) = paste0("b_nino34", colnames(basis_Nino34SSTA))
colnames(basis_EVI) = paste0("b_EVI", colnames(basis_EVI))
colnames(basis_NDVI) = paste0("b_NDVI", colnames(basis_NDVI))
colnames(basis_NDWI) = paste0("b_NDWI", colnames(basis_NDWI))
colnames(basis_LSD) = paste0("b_LSD", colnames(basis_LSD))
colnames(basis_LSN) = paste0("b_LSN", colnames(basis_LSN))


#######################################



data$nYear<-data$Year-2001
data$nCanton=as.numeric(factor(data$Canton))

df <- data %>% 
  dplyr::select(Canton,nCanton,Year,time,nYear,Month,RR,Cases,Poblacion,OFF,
                cluster6_mvc,cluster6_wl,cluster7_wl) %>%
  mutate(Y = Cases, cluster6_mvc=as.numeric(cluster6_mvc), cluster6_wl=as.numeric(cluster6_wl),
         cluster7_wl=as.numeric(cluster7_wl)) %>%
  mutate(date= lubridate::ymd(paste(Year, Month, "01", sep="-")))

df1<-df
df1$Cases[df1$nYear==20]<-NA


#load distance matrices

load(file = './Data/W.canton.vecinos.RData')
load(file = './Data/W.canton.distancia.tiempo.RData')

# W.Distancias<-W.Distancias/max(W.Distancias)
# W.Tiempo<-W.Tiempo/max(W.Tiempo)


quantile.distance <- quantile(W.Distancias[lower.tri(W.Distancias)],p=0.5)
W.adj1<-W.Distancias
W.adj1[W.Distancias<quantile.distance] <- 1
W.adj1[W.Distancias>quantile.distance] <- 0
diag(W.adj1) <- 0


W.adj2<-1/W.Distancias
diag(W.adj2)<-0
W.adj2<- W.adj2/max(W.adj2)


# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# output plot 
indices_cantones <- c(5,6,7,14,19,23,25,26,27,31)
indices_cantones <- c(18,7,28,1,25,19,32,24,16)


