# Packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(tidyr)
library(xts)


hojas <- excel_sheets(path = './Data/Casos_Dengue_1993-2021.xlsx')
annos_sin_53 <- (1994:2021)[-c(4, 10, 15, 21)]
annos_con_53 <- (1994:2021)[c(4, 10, 15, 21, 27)]

# Extraccion de casos observados en Costa Rica ----
CR_datos <-
  read_excel('./Data/Casos_Dengue_1993-2021.xlsx',
             range = "B4:AE57",
             col_names = T) %>%
  gather(key = 'Year', value = 'Cases', `1993`:`2021`) %>%
  filter(Year != 1993) %>%
  filter(Year %in% annos_con_53 | SEMANA != 53)%>%
  na.omit()

fechastot <- ymd('1994-1-8') + weeks(1:dim(CR_datos)[1]-1) #!!!
CR_datos <- cbind(CR_datos,fechastot)
CR_datos_xts <- xts(x = CR_datos$Cases,order.by = CR_datos$fechastot)

# Extraccion de casos observados a nivel cantonal ----
cantones <- read_csv('./Data/32_cantones.csv')
cantones_datos <- NULL

for(i in 1:length(cantones$Cantones)) {
  show(i)
  canton <-
    read_excel(
      './Data/Casos_Dengue_1993-2021.xlsx',
      range = "B4:AE57",
      sheet = cantones$Cantones[i],
      col_names = T
    ) %>%
    gather(key = 'Year', value = 'Cases', `1993`:`2021`) %>%
    filter(Year != 1993) %>%
    filter(Year %in% annos_con_53 | SEMANA != 53) %>%
    mutate(Canton = cantones$Cantones[i])%>%
    na.omit()
  
  fechastot <- ymd('1994-1-8') + weeks(1:dim(canton)[1] - 1)
  canton <- cbind(canton, fechastot)
  cantones_datos <-  bind_rows(cantones_datos, canton)
}

cantones_datos$Canton <- as.factor(cantones_datos$Canton)
cantones_datos$Year <- as.numeric(cantones_datos$Year)

# Filtrado de datos cantonales desde el 2000 y quitando las últimas semanas
# del año en curso o semanas faltantes
cantones_datos <- cantones_datos %>%
  filter(Year >= 2000) %>%
  group_by(Canton) %>%
  dplyr::mutate(semanatot = 1:n()) %>%
  ungroup()

# Cambio de datoscantones a frecuencia mensual
cantones_datos <- cantones_datos %>%
  dplyr::mutate(Month = month(fechastot)) %>%
  dplyr::select(Year, Cases, Canton, Month) %>%
  group_by(Year, Month, Canton) %>%
  summarise(Cases = sum(Cases)) %>%
  ungroup() %>%
  mutate(Year = as.numeric(Year)) %>%
  arrange(Canton,Year,Month)

prop0s <- cantones_datos %>% group_by(Canton,Cases) %>%
   summarise(Totales=n()) %>% ungroup %>%
  group_by(Canton) %>%
  mutate(GTotal=sum(Totales)) %>%
  filter(Cases==0) %>% mutate(Por0=Totales/GTotal)


# Carga del índice SST ----
# Este índice debe ser actualizado constantemente a través de la página:
datosSST <-
  read_table('https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii')

datosSST <- datosSST %>% dplyr::select(
  Year = YR,
  Month = MON,
  Nino12SSTA = ANOM,
  Nino3SSTA = ANOM_1,
  Nino4SSTA = ANOM_2,
  Nino34SSTA = ANOM_3
) %>%
  filter(Year >= 2000, Year <= 2021)


# Carga del índice TNA ----
# Este índice debe ser actualizado constantemente a través de la página:
datosTNA <-
  read_table(
    'https://www.esrl.noaa.gov/psd/data/correlation/tna.data'
    ,
    skip = 1,
    n_max = 73
  )
colnames(datosTNA) <- c('Year', seq(1, 12))

datosTNA <-
  datosTNA %>% pivot_longer(-Year, names_to = 'Month', values_to = 'TNA') %>%
  filter(Year >= 2000, Year <= 2021) %>% mutate(Month = as.numeric(Month))%>%
  filter(TNA>=-90)


# Unión del número de casos cantonales e índices SST y TNA (se remueven NAs
# en el número de casos):
cantones_datos <- cantones_datos %>%
  inner_join(datosSST) %>%
  inner_join(datosTNA) %>%
  mutate(Cases=replace_na(Cases,0))


# Cambio de frecuencia semanal a mensual en casos de CR
CR_datos <- CR_datos %>% 
  mutate(Year=as.numeric(Year)) %>% 
  filter(Year>=2000) %>%
  mutate(Month=month(fechastot)) %>% 
  group_by(Year,Month) %>%
  summarise(CasesCR=sum(Cases)) %>% ungroup()


# Determinación de las constantes de riesgo relativo mensual por cantón----

# Extracción de unidades de área y centros de distritos a partir del 
# Nomenclator http://sistemas.inec.cr/sitiosen/sitiosen/FrmGeografico.aspx

nomenclatura <-
  read_xlsx(path = 'Data/Nomenclator_U-1.xlsx', range = 'A6:X11964')
nomenclatura <-
  nomenclatura %>% dplyr::select(
    DTA = DTA_II,
    CCanton = `Código Canton`,
    Canton = `Nombre Canton`,
    CDistrito = `Código Distrito`,
    Distrito = `Nombre Distrito`,
    CTipo = `Codigo Tipologia`,
    Tipo = Tipologia,
    Latitud,
    Longitud
  ) %>%
  mutate(
    CDistrito = as.numeric(CDistrito),
    CTipo = as.numeric(CTipo),
    Latitud = as.numeric(str_replace(Latitud, ',', '.')),
    Longitud = as.numeric(str_replace(Longitud, ',', '.'))
  ) %>%
  group_by(DTA) %>% arrange(CTipo) %>% slice(1) %>% ungroup()

etiquetas <- nomenclatura %>% dplyr::select(-CTipo,-Tipo,-Latitud,-Longitud)
#write_csv(etiquetas,path = './Data/etiquetas.csv')

# El archivo 'Poblacion_Distrital.xlsx' se crea al unir las etiquetas con las 
# proyecciones de población distrital del INEC 
  # http://www.inec.go.cr/poblacion/estimaciones-y-proyecciones-de-poblacion 

poblacion_distrital <- read_xlsx('Data/Poblacion_Distrital.xlsx')

poblacion_distrital <- poblacion_distrital %>% 
  right_join(cantones,by = 'CCanton') %>% 
  pivot_longer(starts_with('20'),names_to = 'Year',values_to = 'Poblacion') %>%
  dplyr::select(CCanton,Canton,Year,Poblacion,DTA) 

# Cálculo de la población cantonal a nivel anual
poblacion_cantonal <- poblacion_distrital %>%
  group_by(CCanton,Year) %>%
  summarise(Poblacion=sum(Poblacion)) %>% ungroup() %>%
  mutate(Month=1) %>% mutate(Year=as.numeric(Year))


# Cálculo de la población cantonal a nivel mensual (interpolación lineal) 
fechas_base <- CR_datos %>% dplyr::select(Year,Month) %>%
  bind_rows(data.frame(Year=c(rep(2021,7),2022),Month=c(6:12,1))) %>% 
  expand_grid(CCanton=cantones$CCanton) 

poblacion_cantonal <- poblacion_cantonal %>%
  right_join(fechas_base) %>%
  arrange(CCanton,Year,Month) %>%
  mutate(Poblacionap = na.approx(Poblacion,rule=2))%>%
  dplyr::select(Year,Month,CCanton,Poblacion=Poblacionap)

fechas_base_CR <- fechas_base %>% filter(CCanton==101) %>%
  mutate(CCanton=1)

poblacion_CR <- read_xlsx('Data/Poblacion_Nacional.xlsx') %>%
  mutate(Month=1,CCanton=1) %>% 
  right_join(fechas_base_CR) %>%
  arrange(Year,Month) %>%
  mutate(Poblacionap = na.approx(PoblacionCR,rule=2))%>%
  dplyr::select(Year,Month,PoblacionCR=Poblacionap)

# Cálculo de la constante del riesgo relativo por cantón, mes, año  
constante_RR <- poblacion_cantonal %>%
  left_join(poblacion_CR ,by = c('Year','Month')) %>% 
  left_join(CR_datos,by = c('Year','Month')) %>% filter(!is.na(CasesCR)) %>%
  mutate(constRR = (1/Poblacion)/(CasesCR/PoblacionCR))

# Cálculo de la población distrital a nivel mensual (interpolación lineal) 
fechas_base_dist <- CR_datos %>% dplyr::select(Year,Month) %>%
  bind_rows(data.frame(Year=c(rep(2021,7),2022),Month=c(6:12,1))) %>% 
  expand_grid(DTA=unique(poblacion_distrital$DTA))  

poblacion_distrital <- poblacion_distrital %>%
  mutate(Year=as.numeric(Year)) %>% mutate(Month=1) %>% 
  dplyr::select(Year,Poblacion,DTA,Month)

poblacion_distrital <- poblacion_distrital %>%
  right_join(fechas_base_dist,by = c("Year", "DTA", "Month")) %>%
  arrange(DTA,Year,Month) %>%
  mutate(Poblacionap = na.approx(Poblacion,rule=2))%>%
  dplyr::select(Year,Month,DTA,Poblacion=Poblacionap)

# Tabla localizaciones de distritos (MODIS)

localizaciones_MODIS <- nomenclatura %>%
  filter(CCanton %in% cantones$CCanton) %>%
  dplyr::select(CCanton,DTA,Latitud,Longitud)

#save(localizaciones_MODIS,file = './Data/localizaciones_MODIS.RData')

# Carga y tratamiento de datos de MODIS ----

load('./Data/MODISData.RData')

# Extracción de información y unión de bases
baseNDVI <- baseNDVI %>% mutate(variable='NDVI')
baseEVI <- baseEVI %>% mutate(variable='EVI')
baseMIR <- baseMIR %>% mutate(variable='MIR')
baseNIR <- baseNIR %>% mutate(variable='NIR')
baseLSTDia <- baseLST_Dia %>% mutate(variable='LSD')
baseLSTNoche <- baseLST_Noche %>% mutate(variable='LSN')

baseMODIS <-
  bind_rows(baseNDVI, baseEVI, baseMIR, baseNIR, baseLSTDia, baseLSTNoche)
rm(baseNDVI, baseEVI, baseMIR, baseNIR, baseLSTDia, baseLSTNoche)

baseMODIS <- baseMODIS %>% 
  dplyr::select(scale,DTA=site,calendar_date,value,variable)

baseMODIS <- baseMODIS %>%
  mutate(scale = as.numeric(scale), value = as.numeric(value)) %>%
  mutate(ValorMODIS = scale * value,fecha=as_date(calendar_date)) %>%
  mutate(
      ano = year(fecha),
      Month = month(fecha)
    ) %>% dplyr::select(Year=ano,Month,ValorMODIS,variable,DTA) %>% 
  group_by(Year,Month,variable,DTA) %>%
  mutate(ValorMODIS=mean(ValorMODIS)) %>% ungroup() %>% distinct() %>%
  mutate(CCanton=str_extract(DTA, "^.{3}")) %>%
  mutate(CCanton=as.numeric(CCanton),DTA=as.numeric(DTA))

# Base con información MODIS mensual, por distrito con variables
# por columnas, y rellenado de datos con interpolación lineal (pocos casos)
baseMODIS_wide <- baseMODIS %>%
  pivot_wider(names_from = variable,values_from = ValorMODIS) %>%
  mutate(NDWI=(NIR-MIR)/(NIR+MIR)) %>%
  arrange(DTA,Year,Month) %>%
  mutate(NDVI = na.approx(NDVI,rule=2),NDWI = na.approx(NDWI,rule=2),
         EVI = na.approx(EVI,rule=2)) 
  

# Conversión de variables MODIS a frecuencia cantonal, usando
# la población relativa distrital como ponderador

poblacion_distrital <- poblacion_distrital %>% 
  mutate(CCanton=as.numeric(str_extract(DTA, "^.{3}"))) %>%
  group_by(CCanton,Year,Month) %>% distinct() %>%
  mutate(Poblacionrel=Poblacion/sum(Poblacion)) %>% ungroup()  %>%
  dplyr::select(-CCanton)

base_MODIS_cant <- baseMODIS_wide %>% 
  right_join(poblacion_distrital,by=c('Year','Month','DTA')) %>%
  mutate(EVIrel=EVI*Poblacionrel,NDVIrel=NDVI*Poblacionrel,NDWIrel=NDWI*Poblacionrel,LSDrel=LSD*Poblacionrel,LSNrel=LSN*Poblacionrel) %>%
  group_by(Year,Month,CCanton) %>% summarise(EVI=sum(EVIrel),NDVI=sum(NDVIrel),NDWI=sum(NDWIrel),LSD=sum(LSDrel),LSN=sum(LSNrel))%>%
  arrange(CCanton,Year,Month) %>% ungroup() %>%
  drop_na()


# Unión de casos, MODIS, índices y constante de Riesgo Relativo ----
cantones <- cantones %>% dplyr::select(Canton=Cantones,CCanton)
cantones_datos <- cantones_datos %>%
  left_join(cantones,by='Canton') %>%
  filter(Year>2000 | Month >1)

datos_totales <- cantones_datos %>% 
  left_join(constante_RR,by=c('Year','Month','CCanton')) %>%
  left_join(base_MODIS_cant,by=c('Year','Month','CCanton')) %>%
  mutate(RR = Cases*constRR,OFF=log(1/constRR)) #%>%
  #drop_na()


# Carga y tratamiento de datos de CHIRPS ----
load('Data/datos_precipitacion.RData')

localizaciones_MODIS <- localizaciones_MODIS %>%
  mutate(id=1:n())

datos_precipitacion <- datos_precipitacion %>%
  left_join(localizaciones_MODIS,by = 'id') %>%
  group_by(CCanton,date) %>% summarise(Precip=sum(chirps)) %>%
  rename(Fecha=date) %>% ungroup() %>%
  mutate(Year=year(Fecha),Month=month(Fecha)) %>%
  group_by(Year,Month,CCanton) %>%
  summarise(Precip_t=sum(Precip)) %>%
  ungroup()

datos_totales <- datos_totales %>% 
  left_join(datos_precipitacion,by = c('CCanton','Year','Month'))


save(datos_totales,file = './Data/datos_totales.RData')
