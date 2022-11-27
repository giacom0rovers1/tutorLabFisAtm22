###IMPORTARE I DATI
# Clear the workspace
rm(list=ls())
# carico pacchetti
library(tidyverse)
library(lubridate)
library(zoo)
library(terra)
library(rgdal)
library(patchwork)
library(ggsci)
library(rasterVis)


# definizione dei nomi dei files e delle cartelle
rootfolder <- "C:/Users/aless/Desktop/esp_3/"
datafolder <- "C:/Users/aless/Desktop/esp_3/dati/"
HSAF.folder <- paste0(datafolder, "SAT/")
RG.folder   <- paste0(datafolder, "PLV/")
GPM.folder   <- paste0(datafolder, "GPM/")
ERA5.folder   <- paste0(datafolder, "ERA/")
figfolder  <- "figure/"
resfolder  <- "risultati/"

# intervallo temporale di interesse
sel_START <- "201610200000"
sel_END   <- "201610220000"
#trasformo le due stringhe in variabili temporali, poi costruisco la sequenza di tutte le ore tra inizio e fine, ritrasformo in stringhe. Ottengo un array di riferimento della serie temporale
sel.times <- format(seq(ymd_hm(sel_START)+hours(1), 
                        ymd_hm(sel_END), 
                        by="1 hour"), "%Y%m%d%H%M")
# fattori di conversione delle unità di misura (alcuni prodotti sono da convertire in mm/h)
CorrectionFactor_HSAF = 3600 # kg/m2*s --> mm/h
CorrectionFactor_ERA5 = 1000 # m/h --> mm/h
###SOGLIA
THR <- 0.1 # mm/h
#importo shapefile (vettoriali di confini geografici) che servono per avere: confini da plottare come riferimento grafico;  un dominio comune a tutti i sensori con cui mascherare tutti i raster (cioè ritagliare lungo i bordi)
GlobalBG <- vect(paste0(datafolder,"World_Countries_(Generalized)"))
#seleziono il poligono dell'Italia
ItalyBG  <- GlobalBG[112]
#costruisco griglia comune "MyGrid" per avere prodotti con stessi confini e riferiti agli stessi punti ( i prodotti di partenza hanno localizzazioni e risoluzioni diverse). Imposto manualmente le proprietà della griglia.
# Select properties of the common grid ("MyGrid"):
MyGrid_res = 0.1
MyGrid_extLON = c( 5, 20)
MyGrid_extLAT = c(35, 50)
MyGrid <- rast(ncol = round(diff(MyGrid_extLON)/MyGrid_res),
               xmin = MyGrid_extLON[1], 
               xmax = MyGrid_extLON[2], 
               nrow = round(diff(MyGrid_extLAT)/MyGrid_res), 
               ymin = MyGrid_extLAT[1], 
               ymax = MyGrid_extLAT[2])
#bounding box: creo un poligono che include per intero l'italia e il suo intorno per praticità di visualizzazione.
#creo prima un dataframe con le info necessarie e poi un poligono
bbox.lonlat <- cbind(
  id = 1, 
  part = 1, 
  lon = c(MyGrid_extLON, rev(MyGrid_extLON)), 
  lat = c(rep(MyGrid_extLAT[1], 2), rep(MyGrid_extLAT[2], 2))
)
bbox.polyg <- vect(bbox.lonlat, type="polygons")

#####LETTURA RASTER DA FILE#####
ERA5 <- rast(paste0(resfolder, "ERA5_20_21_10.tif"))
RG   <- rast(paste0(resfolder, "RG_20_21_10.tif"))
GPM <- rast(paste0(resfolder, "GPM_20_21_10.tif"))
HSAF <- rast(paste0(resfolder, "HSAF_20_21_10.tif"))
names(RG)[which(! names(RG) %in% names(HSAF))]

new.times <- names(RG)[which(names(RG) %in% names(HSAF))]
RG   <-   RG[[new.times]]
ERA5 <- ERA5[[new.times]] 
GPM  <-  GPM[[new.times]]

#####CONFINI REGIONE SICILIA#####
Regioni <- vect(paste0(datafolder,"Reg01012022_g"))
#seleziono il poligono dell'Italia
Sicilia  <- Regioni[19]
plot(Sicilia)

##applico maschera per considerare dati solo sulla Regione Sicilia
ERA5.mask <- mask(ERA5 >= THR, Sicilia)
GPM.mask  <- mask(GPM  >= THR, Sicilia)
RG.mask   <- mask(RG   >= THR, Sicilia)
HSAF.mask <- mask(!is.na(HSAF), Sicilia)

#####FUNZIONI PER STATISTICA#####
statsFun <- function(Var, Ref){
  
  avgR <- mean(Ref, na.rm=T)
  SMP  <- length(na.omit(Var))
  COR  <- cor(Var, Ref, use = "pairwise.complete.obs")
  ME   <- mean(Var - Ref, na.rm = T)
  MAE  <- mean(abs(Var - Ref), na.rm = T)
  RMSE <- sqrt(mean((Var - Ref)^2, na.rm = T))
  CV   <- RMSE/avgR
  nME  <- ME/avgR
  nMAE <- MAE/avgR
  MAX  <- max(Var, na.rm = T)
  
  table <- data.frame(SMP, avgR, nME, nMAE, CV, COR, ME, MAE, RMSE, MAX)
  
  return(round(table, 2))
  
}

statsFun_CAT <- function(Var, Ref){
  cat = data.frame(Var = Var==1, 
                   Ref = Ref==1)
  
  xtab <- xtabs(~ Ref + Var, data = cat, na.action = "na.omit")
  
  hits <- xtab['TRUE' ,'TRUE' ]  #[2,2]#
  crej <- xtab['FALSE','FALSE']  #[1,1]#
  miss <- xtab['TRUE' ,'FALSE']  #[2,1]#
  fala <- xtab['FALSE','TRUE' ]  #[1,2]#
  
  FAR  <- fala/(hits + fala)
  BIAS <- (hits + fala)/(hits + miss)
  POD  <-  hits/(hits + miss)       
  CSI  <-  hits/(hits + miss + fala)  
  SMP  <- sum(xtab)
  PC   <- (hits+crej)/SMP
  
  ACC  <-  (hits + crej)/(hits + fala + miss + crej)
  POFD <- fala/(fala + crej) 
  
  rand <- (hits+miss)*(hits+fala)/SMP
  ETS  <- (hits - rand)/(hits + miss + fala - rand)
  
  table <- data.frame(SMP, POD, PC, FAR, BIAS, ETS, CSI, ACC, POFD)
  
  return(round(table, 2))
}

#####AGGREGAZIONE DATI#####
data <- tibble(
  RainGauges = as.vector(values(clamp(RG, THR, values=F))),
  ERA5 = as.vector(values(clamp(ERA5, THR, values=F))),
  GPM = as.vector(values(clamp(GPM, THR, values=F))),
  HSAF = as.vector(values(clamp(HSAF, THR, values=F)))
)  %>% na.omit()

data_pivot <- data %>% pivot_longer(
  cols = c("RainGauges", "ERA5", "GPM", "HSAF"),
  names_to = "Prodotti",
  values_to = "Precipitazione (mm/h)")

data_CAT <- tibble(
  RainGauges = as.vector(values(RG.mask)),
  ERA5 = as.vector(values(ERA5.mask)),
  GPM = as.vector(values(GPM.mask)),
  HSAF = as.vector(values(HSAF.mask))
) %>% na.omit()


#####SCATTERPLOTS DI BASE#####

#ERA5 & RAINGAUGES
scatter_ERA5 <- data %>%
  ggplot(aes(RainGauges, ERA5)) +
  geom_bin_2d(binwidth = 0.2) +
  coord_fixed(xlim = c(0, 20),
              ylim = c(0, 20)) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.2) +
  scale_fill_binned(type = "viridis", 
                    trans = "log", 
                    breaks = c(2, 5, 10, 25, 50, 100, 200)) +
  theme_bw()

#HSAF & RAINGAUGES
scatter_HSAF <- data %>%
  ggplot(aes(RainGauges, HSAF)) +
  geom_bin_2d(binwidth = 0.2) +
  coord_fixed(xlim = c(0, 20),
              ylim = c(0, 20)) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.2) +
  scale_fill_binned(type = "viridis", 
                    trans = "log", 
                    breaks = c(2, 5, 10, 25, 50, 100, 200)) +
  theme_bw()

#GPM & RAINGAUGES
scatter_GPM <- data %>%
  ggplot(aes(RainGauges, GPM)) +
  geom_bin_2d(binwidth = 0.2) +
  coord_fixed(xlim = c(0, 20),
              ylim = c(0, 20)) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.2) +
  scale_fill_binned(type = "viridis", 
                    trans = "log", 
                    breaks = c(2, 5, 10, 25, 50, 100, 200)) +
  theme_bw()

#scatter <- scatter_ERA5 + scatter_HSAF + scatter_GPM + plot_annotation(title="Scatterplot")
show(scatter_ERA5)
show(scatter_HSAF)
show(scatter_GPM)
#print(scatter)
#ggsave("figure/coperture.png" , coperture)


#####STATISTICHE DI CONFRONTO - Continue#####
#Confronto ERA5 & RG
stats_cont_ERA5 <- statsFun(data$ERA5, data$RainGauges)
stats_cont_ERA5

#Confronto HSAF & RG
stats_cont_HSAF <- statsFun(data$HSAF, data$RainGauges)
stats_cont_HSAF

#Confronto GPM & RG
stats_cont_GPM <- statsFun(data$GPM, data$RainGauges)
stats_cont_GPM

#####STATISTICHE DI CONFRONTO - Categoriche#####
#Confronto ERA5 & RG
stats_CAT_ERA5 <- statsFun_CAT(data_CAT$ERA5, data_CAT$RainGauges)
stats_CAT_ERA5

#Confronto HSAF & RG
stats_CAT_HSAF <- statsFun_CAT(data_CAT$HSAF, data_CAT$RainGauges)
stats_CAT_HSAF

#Confronto GPM & RG
stats_CAT_GPM <- statsFun_CAT(data_CAT$GPM, data_CAT$RainGauges)
stats_CAT_GPM

#####COMPORTAMENTO INDICI CATEGORICI AL VARIARE DELLA SOGLIA
THR_vec<-c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
#, 13, 14, 15)
table <- data.frame()
#c'è un problema con la soglia 0.1 e da 13 in poi
for (i in THR_vec){
  ERA5.mask_temp <- mask(ERA5 >= i, ItalyBG)
  GPM.mask_temp  <- mask(GPM  >= i, ItalyBG)
  RG.mask_temp  <- mask(RG   >= i, ItalyBG)
  HSAF.mask_temp <- mask(HSAF >= i, ItalyBG)
  
  data_CAT_temp <- tibble(
    RainGauges = as.vector(values(RG.mask_temp)),
    ERA5 = as.vector(values(ERA5.mask_temp)),
    GPM = as.vector(values(GPM.mask_temp)),
    HSAF = as.vector(values(HSAF.mask_temp))
  ) %>% na.omit()
  
  table <- rbind(table,cbind(Soglia=c(i), Prodotto=c(factor("ERA5")), statsFun_CAT(data_CAT_temp$ERA5, data_CAT_temp$RainGauges)))
  table <- rbind(table,cbind(Soglia=c(i), Prodotto=c(factor("HSAF")), statsFun_CAT(data_CAT_temp$HSAF, data_CAT_temp$RainGauges)))
  table <- rbind(table,cbind(Soglia=c(i), Prodotto=c(factor("GPM")), statsFun_CAT(data_CAT_temp$GPM, data_CAT_temp$RainGauges)))
  
}

graf_POD <- table %>% ggplot(aes(Soglia, POD, group = Prodotto, colour=Prodotto, shape=Prodotto)) + 
  geom_point() + 
  geom_line()+
  labs(x="Soglia dry-wet (mm/h)")+
  theme_bw()
graf_FAR <- table %>% ggplot(aes(Soglia, FAR, group = Prodotto, colour=Prodotto, shape=Prodotto)) + 
  geom_point() + 
  geom_line()+
  labs(x="Soglia dry-wet (mm/h)")+
  theme_bw()
graf_CSI <- table %>% ggplot(aes(Soglia, CSI, group = Prodotto, colour=Prodotto, shape=Prodotto)) + 
  geom_point() + 
  geom_line()+
  labs(x="Soglia dry-wet (mm/h)")+
  theme_bw()
graf_BIAS <- table %>% ggplot(aes(Soglia, BIAS, group = Prodotto, colour=Prodotto, shape=Prodotto)) + 
  geom_point() + 
  geom_line()+
  labs(x="Soglia dry-wet (mm/h)")+
  theme_bw()
graf_PC <- table %>% ggplot(aes(Soglia, PC, group = Prodotto, colour=Prodotto, shape=Prodotto)) + 
  geom_point() + 
  geom_line()+
  labs(x="Soglia dry-wet (mm/h)")+
  theme_bw()
graf_ETS <- table %>% ggplot(aes(Soglia, ETS, group = Prodotto, colour=Prodotto, shape=Prodotto)) + 
  geom_point() + 
  geom_line()+
  labs(x="Soglia dry-wet (mm/h)")+
  theme_bw()

graf_INDICI<-graf_POD+graf_FAR+graf_CSI+graf_ETS+graf_BIAS+graf_PC+guide_area()+plot_layout(guides="collect")

show(graf_INDICI)


#####PRODURRE MAPPE SIMULTANEE DI PRODOTTI DIVERSI#####
plotMultiMap <- function(HSAF, GPM, ERA5, RG, layerNames){
  
  nCols <- 2
  nRows <- 2
  par(mfrow=c(nRows, nCols))
  
  # list with details for legend
  plg = list(
    title = expression(bold("Average rain\n rate (mm/h)")),
    title.cex = 0.7,
    cex = 0.7,
    shrink=0
  )
  
  # margins: bottom, left, top, and right
  mar = c(4.5, 4.5, 2.5, 3.8)
  
  for(lname in layerNames){
    for(rname in c("HSAF", "GPM", "ERA5", "RG")){
      raster <- get(rname)
      rastLayer <- raster[[lname]]
      
      plot(
        clamp(rastLayer, 0.1, values=F), 
        main = paste(rname, " - ", ymd_hm(lname)),
        col  = rev(rainbow(295)[1:200]),  
        xlab = "Longitude (°E)",
        ylab = "Latitude (°N)",
        xlim = c(12, 16),
        ylim = c(36, 39),
        range=c(0, 17),
        #        plg  = plg,
        mar = mar,
        panel.first = {
          plot(
            crop(GlobalBG, MyGrid), 
            col="gray90", 
            add=T
          )
          
          grid()
        }
      )
    }
  }
}
range_temp <- format(seq(ymd_hm(sel_START)+hours(1), 
                         ymd_hm(sel_END), 
                         by="1 hour"), "%Y%m%d%H%M")

for (i in range_temp){
  plotMultiMap(HSAF, GPM, ERA5, RG, i)
}
