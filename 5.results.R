library(tidyverse)
library(kableExtra)

metrica.canton.in <- read.csv("results/metricas_Canton_INLA_in.csv")
metrica.canton.out <- read.csv("results/metricas_Canton_INLA_out.csv")

metrica.DTW6.in <- read.csv("results/metricas_cluster_DTW6_INLA_in.csv")
metrica.DTW6.out <- read.csv("results/metricas_cluster_DTW6_INLA_out.csv")

metrica.WL6.in <- read.csv("results/metricas_cluster_wl6_INLA_in.csv")
metrica.WL6.out <- read.csv("results/metricas_cluster_wl6_INLA_out.csv")


tabla.in <- metrica.canton.in %>% left_join(metrica.DTW6.in,by="Canton") %>% 
      left_join(metrica.DTW6.in,by="Canton")

tabla.out <- metrica.canton.out %>% left_join(metrica.DTW6.out,by="Canton") %>% 
  left_join(metrica.DTW6.out,by="Canton")
