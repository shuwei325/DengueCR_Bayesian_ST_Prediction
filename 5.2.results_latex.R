library(tidyverse)
library(kableExtra)
library(data.table)


# metrics -----------------------------------------------------------------

# lineal

fitting <- read.csv("results_out3mes/table_fitting.csv")

fitting1 <- fitting[c(1,14,15,16,18,19,20),] %>% mutate(spatial=recode(Model,
                                   'formula.iid'='Independent',
                                   'formula.1.1a'='CAR',
                                   'formula.1.2a'='proper CAR',
                                   'formula.1.3a'='BYM',
                                   'formula.2.1a'='CAR',
                                   'formula.2.2a'='proper CAR',
                                   'formula.2.3a'='BYM'),
                                   distance=c("",rep(c("Neighbor","Distance"),each=3)),
                                   DLNM = rep("Linear",7))

fitting2 <- fitting1 %>% select(DLNM,distance,spatial,DIC, logscore)

fitting2 %>% kbl(digits=4,row.names=FALSE)%>%
  kable_minimal() %>%
  kable_paper(full_width = F)

fitting2 %>% kbl(digits=4,format="latex",row.names=FALSE)%>%
  kable_minimal() %>%
  kable_paper(full_width = F)


# non-lineal

fitting.nonlinear <- read.csv("results_nolineal/table_fitting.csv")
fitting.nonlinear1 <- fitting.nonlinear[c(1,14,15,16,18,19,20),] %>% mutate(spatial=recode(Model,
                                                                       'formula.iid'='Independent',
                                                                       'formula.1.1a'='CAR',
                                                                       'formula.1.2a'='proper CAR',
                                                                       'formula.1.3a'='BYM',
                                                                       'formula.2.1a'='CAR',
                                                                       'formula.2.2a'='proper CAR',
                                                                       'formula.2.3a'='BYM'),
                                                        distance=c("",rep(c("Neighbor","Distance"),each=3)),
                                                        DLNM = rep("Non-linear",7))

fitting.nonlinear2 <- fitting.nonlinear1 %>% select(DLNM,distance,spatial,DIC, logscore)

fitting.all <- fitting2 %>% bind_rows(fitting.nonlinear2)



fitting.all %>% kbl(digits=4,row.names=FALSE)%>%
  kable_minimal() %>%
  kable_paper(full_width = F)

fitting.all %>% kbl(digits=4,format="latex",row.names=FALSE)%>%
  kable_minimal() %>%
  kable_paper(full_width = F)

fitting.all %>% arrange(DIC) %>% kbl(digits=4,row.names=FALSE)%>%
  kable_minimal() %>%
  kable_paper(full_width = F)

fitting.all %>% arrange(logscore) %>% kbl(digits=4,row.names=FALSE)%>%
  kable_minimal() %>%
  kable_paper(full_width = F)


#dd


metrica.iid.in <- read.csv("results_out3mes/formula.iid.metricas.in.csv")
metrica.iid.out <- read.csv("results_out3mes/formula.iid.metricas.out.csv")

metrica.1.2a.in <- read.csv("results_out3mes/formula.1.2a.metricas.in.csv")
metrica.1.2a.out <- read.csv("results_out3mes/formula.1.2a.metricas.out.csv")


tabla.in <- metrica.iid.in %>% left_join(metrica.1.2a.in,by="Canton")

tabla.out <- metrica.iid.out %>% left_join(metrica.1.2a.out,by="Canton")

options(scipen=999)


tabla.in %>% kbl(digits=2)%>%
  kable_minimal() %>%
  kable_paper(full_width = F) %>%
  add_header_above(c(" " = 1, "iid" = 2, 
                     "spatial" = 2))

tabla.in %>% kbl(digits=4,format="latex",row.names=FALSE)%>%
  kable_minimal() %>% 
  kable_paper(full_width = F) %>%
  add_header_above(c(" " = 1, "iid" = 2, 
                     "spatial" = 2))


tabla.out %>% kbl(digits=2)%>%
  kable_minimal() %>%
  kable_paper(full_width = F) %>%
  add_header_above(c(" " = 1, "iid" = 2, 
                     "spatial" = 2))


tabla.out %>% kbl(digits=4,format="latex",row.names=FALSE)%>%
  kable_minimal() %>% 
  kable_paper(full_width = F) %>%
  add_header_above(c(" " = 1, "iid" = 2, 
                     "spatial" = 2))


