library(dplyr)
library(pscl)
library(VGAM)
library(ggplot2)
library(tidyverse)
library(GGally)
library(gamlss)
library(astsa)
library(psych)
library(lubridate)
library(tidyverse)
library(reshape2)
library(dtw)
library(gridExtra)
library(dtwclust)

# 0. reading data ------------------------------------------------------------

load(file = './Data/datos_totales.RData')
str(datos_totales)
dim(datos_totales)
table(datos_totales$Year,datos_totales$Month)


# 1. Preparing data -------------------------------------------------------

#Filter the data from 2000-02 to 2021-03
datos_totales <- datos_totales %>% 
                      dplyr::filter(Year!=2021 | Month %in% c(1,2,3,3)) 
                        
# creating time index
datos_totales <-datos_totales %>%
                      mutate(time=make_date(Year,Month,15))

#datos_totales1<-na.omit(datos_totales)

dim(datos_totales)
names(datos_totales)
colSums(is.na(datos_totales))

# creating logCases, logRR and logPrecip
datos_totales2 <- datos_totales %>% 
  dplyr::mutate(logCases=log(Cases+0.5),logRR=log(RR+0.5),logPrecip=log(Precip_t+0.5))%>% 
  dplyr::select(Canton,time,Cases,logCases,RR,logRR,Precip_t,logPrecip,Nino12SSTA,Nino3SSTA,Nino4SSTA,Nino34SSTA,TNA,
                Poblacion, PoblacionCR,
                EVI,NDVI,NDWI,LSD,LSN,OFF) 

names(datos_totales2)
#Nino
pairs.panels(datos_totales2[,c(6,9:13)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
pairs.panels(datos_totales2[,c(6,16:20)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

# 2. Converting data.frame to wide format ---------------------------------


# 2.1. logCases -----------------------------------------------------------

logCases_wide <- dcast(datos_totales2, Canton ~ time, value.var="logCases")
logCases_wide1 <- logCases_wide %>% dplyr::select(-Canton)
dim(logCases_wide1)
distMatrix_logCases <- dist(logCases_wide1, method="DTW")

logCases_hc <- hclust(distMatrix_logCases, method="average")

plot(logCases_hc,  main="")
rect.hclust(logCases_hc, k = 6, border = 5)


cluster6_logCases <- cutree(logCases_hc, k = 6)
logCases_wide$cluster5_logCases <- as.factor(cluster5_logCases)
logCases_wide$cluster6_logCases <- as.factor(cluster6_logCases)

results_logCases <- logCases_wide %>% dplyr::select(Canton,cluster5_logCases, cluster6_logCases)


# 2.2. logRR --------------------------------------------------------------
logRR_wide <- dcast(datos_totales2, Canton ~ time, value.var="logRR")
logRR_wide1 <- logRR_wide %>% dplyr::select(-Canton)
dim(logRR_wide1)

distMatrix_logRR <- dist(logRR_wide1, method="DTW")

logRR_hc <- hclust(distMatrix_logRR, method="average")

plot(logRR_hc,  main="")
rect.hclust(logRR_hc, k = 5, border = 4) # 4 clusters
rect.hclust(logRR_hc, k = 6, border = 5)

cluster5_logRR <- cutree(logRR_hc, k = 5)
cluster6_logRR <- cutree(logRR_hc, k = 6)
logRR_wide$cluster5_logRR <- as.factor(cluster5_logRR)
logRR_wide$cluster6_logRR <- as.factor(cluster6_logRR)

results_logRR <- logRR_wide %>% dplyr::select(Canton, cluster5_logRR, cluster6_logRR)

# 2.3. logPrecip ----------------------------------------------------------

logPrecip_wide <- dcast(datos_totales2, Canton ~ time, value.var="logPrecip")
logPrecip_wide1 <- logPrecip_wide %>% dplyr::select(-Canton)
row.names(logPrecip_wide1)<-logPrecip_wide$Canton
dim(logPrecip_wide1)

distMatrix_logPrecip <- dist(logPrecip_wide1, method="DTW")

logPrecip_hc <- hclust(distMatrix_logPrecip, method="average")

plot(logPrecip_hc,  main="")
rect.hclust(logPrecip_hc, k = 5, border = 4)
rect.hclust(logPrecip_hc, k = 6, border = 5)

cluster5_logPrecip <- cutree(logPrecip_hc, k = 5)
cluster6_logPrecip <- cutree(logPrecip_hc, k = 6)
logPrecip_wide$cluster5_logPrecip <- as.factor(cluster5_logPrecip)
logPrecip_wide$cluster6_logPrecip <- as.factor(cluster6_logPrecip)

results_logPrecip <- logPrecip_wide %>% dplyr::select(Canton, cluster5_logPrecip, cluster6_logPrecip)
  
#######################

#Comparación de clustering usando logCases y logRR.
table(results_logCases$cluster6_logCases,results_logRR$cluster6_logRR)
table(results_logPrecip$cluster6_logPrecip,results_logRR$cluster6_logRR)


datos_totales3 <- datos_totales2 %>% 
                          left_join(results_logCases, by= "Canton") %>% 
                          left_join(results_logRR, by= "Canton") %>%
                          left_join(results_logPrecip, by= "Canton")



ggplot(datos_totales3, aes(x = time, y = logCases)) + 
  geom_line(aes(color = cluster5_logCases), size = 1) +
  theme_minimal()

ggplot(datos_totales3, aes(x = time, y = logRR)) + 
  geom_line(aes(color = cluster5_logRR), size = 1) +
  theme_minimal()

ggplot(datos_totales3, aes(x = time, y = logRR)) + 
  geom_line(aes(color = cluster5_logPrecip), size = 1) +
  theme_minimal()

ggplot(datos_totales3, aes(x = time, y = logPrecip)) + 
  geom_line(aes(color = cluster6_logPrecip), size = 1) +
  theme_minimal()

#Group_by

grupos_logRR <- results_logRR %>% dplyr::select(Canton, 
                                   cluster5_logRR, cluster6_logRR) %>% 
            group_by(cluster5_logRR) %>% 
  nest()

grupos_logRR$data[[1]]
grupos_logRR$data[[2]]
grupos_logRR$data[[3]]
grupos_logRR$data[[4]]
grupos_logRR$data[[5]]

grupos <- results_logPrecip %>% dplyr::select(Canton, 
                                          cluster5_logPrecip, cluster6_logPrecip) %>% 
  group_by(cluster5_logPrecip) %>% 
  nest()

grupos_logPrecip$data[[1]]
grupos_logPrecip$data[[2]]
grupos_logPrecip$data[[3]]
grupos_logPrecip$data[[4]]
grupos_logPrecip$data[[5]]


save(results_logRR,results_logPrecip ,file="./Data/clustering.Rdata")


# 3. Using dtwclust -------------------------------------------------------

# 3.2. multivariate ----------------------------------------------------------

names(datos_totales2)
series <- datos_totales2 %>% dplyr::select(Canton,logPrecip,
                                Nino12SSTA,Nino3SSTA,Nino4SSTA,Nino34SSTA,
                                TNA,EVI,NDVI,NDWI,LSD,LSN) %>%
              group_by(Canton) %>% nest()

series.mv <- series$data
names(series.mv) <- series$Canton
labels <- as.factor(series$Canton)
length(series.mv)
class(series.mv)

#DTW_basic

mvc5 <- tsclust(series.mv, type = "hierarchical",
                 distance = "dtw2", k = 5L, trace = TRUE,
                 control = hierarchical_control(method = "all"))
class(mvc5)
length(mvc5)

# Cluster validity indices
CIV5<-sapply(mvc5, cvi, b = labels)
unlist(apply(CIV5,1,which.min))
mvc5a<-mvc5[[1]]

names(mvc5a)
mvc5a@distance
plot(mvc5a)
rect.hclust(mvc5a, k = 5, border = 2) 
mvc5a@cluster
mvc5a@clusinfo

results_mvc5a<-mvc5a@cluster


mvc6 <- tsclust(series.mv, type = "hierarchical",
                distance = "dtw2", k = 6L, trace = TRUE,
                control = hierarchical_control(method = "all",
                                               distmat = mvc5[[1L]]@distmat))
CIV6<-sapply(mvc6, cvi, b = labels)
unlist(apply(CIV6,1,which.min))
mvc6a<-mvc6[[1]]

plot(mvc6a)
rect.hclust(mvc6a, k = 6, border = 2) 
mvc6a@cluster
mvc6a@clusinfo

results_mvc6a<-mvc6a@cluster


mvc4 <- tsclust(series.mv, type = "hierarchical",
                distance = "dtw2", k = 4L, trace = TRUE,
                control = hierarchical_control(method = "all",
                                               distmat = mvc5[[1L]]@distmat))
CIV4<-sapply(mvc4, cvi, b = labels)
unlist(apply(CIV4,1,which.min))
mvc4a<-mvc4[[1]]
results_mvc4a<-mvc4a@cluster

mvc3 <- tsclust(series.mv, type = "hierarchical",
                distance = "dtw2", k = 3L, trace = TRUE,
                control = hierarchical_control(method = "all",
                                               distmat = mvc5[[1L]]@distmat))
CIV3<-sapply(mvc3, cvi, b = labels)
unlist(apply(CIV3,1,which.min))
mvc3a<-mvc3[[1]]
results_mvc3a<-mvc3a@cluster

results_mvc <- data.frame(Canton=names(results_mvc3a),cluster3_mvc=results_mvc3a,cluster4_mvc=results_mvc4a,cluster5_mvc=results_mvc5a,cluster6_mvc=results_mvc6a)


save(results_logRR,results_logPrecip,results_mvc ,file="./Data/clustering.Rdata")


