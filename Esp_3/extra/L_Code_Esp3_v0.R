#-------------------------------------------------------------------------------
#ESPERIENZA 3 09/11
#ANALISI DI CAMPI DI PRECIPITAZIONE STIMATI DA SATELLITE
#-------------------------------------------------------------------------------

# Clear the workspace
rm(list=ls())

# Load packages

packages <- c(
  "tidyverse",
  "lubridate",
  "zoo",
  "terra", #L: Pacchetto fatto apposta per gestire dati geo appartenenti alla terra
  "rgdal"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = T)){
    install.packages(pkg)
  } 
  
  library(pkg, character.only = T)
}

tidyverse_logo()

# Define the working directory (in this case, the path is different depending on whether the code is executed on Linux or Windows, but just for my convenience)


# Definizione dei nomi dei files e delle cartelle--------------------------------
# In questo caso datafolder è la cartella dentro cui si trovano le sottocartelle 
 # resfolder conterrà tutti i dati processati e i risultati, figfolder le figure.
 

# Define the project directories (relative path)
     
datafolder <- "dati3/"
if(!dir.exists(datafolder)){ dir.create(datafolder) }
HSAF.folder <- paste0(datafolder, "SAT/")
RG.folder   <- paste0(datafolder, "PLV/")
GPM.folder   <- paste0(datafolder, "GPM/")
ERA5.folder   <- paste0(datafolder, "ERA/")

figfolder  <- "figure/"

# resfolder  <- "risultati/"

# Create the project directories, if missing
if(!dir.exists(figfolder)){  dir.create(figfolder) }
# if(!dir.exists(resfolder)){  dir.create(resfolder) }

# https://hub.arcgis.com/datasets/2b93b06dc0dc4e809d3c8db5cb96ba69_0/

#L: Avremo a che fare con delle mappe, quindi scarichiamo uno shape file 
GlobalBG <- vect(paste0(datafolder,"World_Countries_(Generalized)")) 
#L: "vect()" li crea già come entità spaziali (latitudine, longitudine e concetti di forma:
 # il software sa già che ha un perimetro, quindi che ha un dentro e un fuori).
ItalyBG  <- GlobalBG[112] #L:112 corrisponde all'Italia

plot(ItalyBG, col="grey")
#Grafico
#Fornisce sia i confini da plottare, sia un dominio con cui mascherare i raster

#GRIGLIA COMUNE ("MY GRID")-----------------------------------------------------
#L: Fare statistica e confronti ha senso solo se ci mettiamo su un campo comune. 
 # Perciò andiamo a costruire una griglia uguale per tutti

#Creo un oggetto SpatRaster sull’Italia

MyGrid <- rast(ncol=250, nrow=150, 
               xmin=0, xmax=25, 
               ymin=35, ymax=50)

MyGrid          #Ho una griglia su tutta l'italia quindi. Quindi se sapessi com'è 
                #disposta potrei selezionare una zona a piacere

#Situazione e obiettivo
#Input: 4 risoluzioni spaziali e temporali diverse Output: stessa griglia, cumulate orarie (o mm/h medi orari)

#Import data--------------------------------------------------------------------
#-------------------------------------------------------------------------------

##ERA5---------------------------------------------------------------------------
CorrectionFactor_ERA5 = 1000 #Fattore correttivo

#L:Andiamo a prendere il nome del file
ERA5.files <- dir(ERA5.folder, pattern=".grib")
length(ERA5.files)

### Import data 
ERA5.raster <- rast(paste0(ERA5.folder, ERA5.files)) * CorrectionFactor_ERA5 #L: Fattore moltiplicativo correttivo
#L: Un raster è già predisposto come una griglia e dei valori associati ai pixel 
 # della griglia ("rast"). Ci pensa lui a gestire e convertire i pixel sulla sup 
 # sferica
dim(ERA5.raster)
# L:Otteniamo un oggetto solo che ha già capito che deve essere stratificato, 
 #dove ogni strato rappresenta un istante temporale

### Assign Date 
ERA5.times_START <- substr(ERA5.files,6,15)
ERA5.times_END   <- substr(ERA5.files,17,26)


ERA5.times <- format(seq(ymd_h(ERA5.times_START), ymd_h(ERA5.times_END), by="1 hour"), "%Y%m%d%H%M")

names(ERA5.raster) <- ERA5.times

ERA5.raster
    #TABELLA
#L: Dando questi estremi (vedi tabella) ottengo una risoluzione di 0.1°

#L: Plot fa cose diverse a seconda di cosa gli arriva
plot(ERA5.raster[["201610012100"]])
plot(ItalyBG, add =T)

    # Grafico
#L: Nel grafico il grigio rappresenta lo 0, mentre il bianco dove non si ha dati

ERA5.MyGrid <- terra::project(ERA5.raster, MyGrid)
#L: Non sto sovrascrivendo, ma li sto proiettando. Project sa fare le interpolazioni
 # e le fa in automatico in lineare, ma con uno specifico comando posso fare fit
 #anche diversi


plot(ERA5.MyGrid[["201610012100"]])
plot(ItalyBG, add =T)


plot(
  #L: "clamp" taglia il raster da 0.2 in giù (è l'analogo di filter?)
  clamp(ERA5.MyGrid[["201610012100"]], 0.2, values=F), 
  main = ymd_hm(names(ERA5.MyGrid[["201610012100"]])),
  col=c("blue","green","yellow","orange","red"),  
  xlab   = "Longitude (°E)",
  ylab   = "Latitude (°N)",
  panel.first = {
    plot(crop(GlobalBG, MyGrid), col="gray90", add=T)
    grid()
  }
)
            #Forse serve installare il pacchetto "rgdal"
      
##H-SAF--------------------------------------------------------------------------
CorrectionFactor_HSAF = 3600 # 1
#importare i dati HSAF e rasterizzarli
### Search GRIB files 
HSAF.files <- dir(HSAF.folder, pattern=".grb")
#HSAF.files <- HSAF.files[substr(HSAF.files,9,12)=="1001"] 
length(HSAF.files)
  
#### Import data (Precipitation Rate [kg/m2 s])
HSAF.list <- list()
#L: Importo i dati dentro una lista come una matrice di dati. Al contrario dei 
 #dataframes non ha bisogno che le entrate siano tutte uniformi
for(i in 1:length(HSAF.files)) {
  HSAF.list[[i]] <- rgdal::readGDAL(paste0(HSAF.folder,HSAF.files[i]), silent = TRUE)
}

length(HSAF.list)

dim(HSAF.list[[1]])
 

### Import Coordinates (unstructured grid)
#L: Ora ho bisogno delle coord di supporto, mentre la info temporale sta dentro 
 #la lista
HSAF.coordinates <- read.table(file(paste0(datafolder,"SAT_coordinates.dat")),
                               sep = " ", na.string = NA, as.is = TRUE, header = FALSE)
dim(HSAF.coordinates)
 

### Assign Date 
HSAF.times <- paste0(substr(HSAF.files,5,12),substr(HSAF.files,14,17))
HSAF.times [c(1,length(HSAF.times))]
    
### Merge Data 
HSAF.df <- HSAF.coordinates[1710000:1,]
for(i in 1:length(HSAF.times)) HSAF.df <- cbind(HSAF.df, matrix(HSAF.list[[i]][[1]], dim(HSAF.list[[1]])[1], 1))
colnames(HSAF.df) <- c("Lon", "Lat", HSAF.times)

head(HSAF.df)

 
dim(HSAF.df)

#L: Visto che è un file molto grande, vado a filtrarlo
HSAF.df_ITA <- HSAF.df %>% filter(
  Lon >= 0,
  Lon <= 30,
  Lat >= 30,
  Lat <= 60
) 

rownames(HSAF.df_ITA) <- 1:nrow(HSAF.df_ITA)

head(HSAF.df_ITA)

 
#L: "rasterize" punto chiave per tornare ad un raster
HSAF.raster <- rasterize(cbind(HSAF.df_ITA$Lon, HSAF.df_ITA$Lat), MyGrid, values=as.array(HSAF.df_ITA[["201610012057"]])) * CorrectionFactor_HSAF

HSAF.raster <- clamp(HSAF.raster, 0, values=F)

names(HSAF.raster) <- "201610012057"
HSAF.raster

 

plot(HSAF.raster)
plot(GlobalBG, add=T)



plot(ItalyBG, col="grey")
plot(HSAF.raster[["201610012057"]], 
     add=T, 
     col=c("blue","green","yellow","orange","red"))

plot(
  HSAF.raster[["201610012057"]], 
  main = names(HSAF.raster[["201610012057"]]),
  col=c("blue","green","yellow","orange","red"),  
  xlab   = "Longitude (°E)",
  ylab   = "Latitude (°N)",
  panel.first = {
    plot(crop(GlobalBG, MyGrid), col="gray90", add=T)
    # plot(ItalyBG, col="gray80", add=T)
    grid()
  }
)


hist(HSAF.raster, "201610012057")

 
##GPM-IMERG---------------------------------------------------------------------
#L:Ha un altro formato ancora e diventa una matrice
### Search IMR files (ASCII) 
GPM.files <- dir(GPM.folder, pattern=".imr")
GPM.files <- GPM.files[substr(GPM.files,5,8)=="1001"]
length(GPM.files)
  
### Import data (Precipitation Rate [mm/h]) 
GPM.list <- list()
for(i in 1:(length(GPM.files)) ) {
  GPM.list[[i]] <- read.table(file(paste0(GPM.folder,GPM.files[i])), skip = 181, nrows = 180,
                              sep = ",", na.string = NA, as.is = TRUE, header = FALSE)
  GPM.list[[i]]  <- as.matrix(GPM.list[[i]][-1])
}


### Import Coordinates (structured grid) 
LatGPM <- read.table(file(paste0(GPM.folder,GPM.files[1])), skip = 361, nrows = 1,
                     sep = ",", na.string = NA, as.is = TRUE, header = FALSE)
temp   <- LatGPM[[2]]
for(i in 3:length(LatGPM)) temp[i-1] <- LatGPM[[i]]
LatGPM <- temp

LonGPM <- read.table(file(paste0(GPM.folder,GPM.files[1])), skip = 362, nrows = 1,
                     sep = ",", na.string = NA, as.is = TRUE, header = FALSE)
temp   <- LonGPM[[2]]
for(i in 3:length(LonGPM)) temp[i-1] <- LonGPM[[i]]
LonGPM <- temp

c(length(LonGPM),length(LatGPM))
  
GPM.coordinates <- 0
for(i in LonGPM) for(j in LatGPM) GPM.coordinates <- c(GPM.coordinates,i,j)
# Data extraction
GPM.coordinates2 <- t(matrix(GPM.coordinates[-1], 2, length(LonGPM)*length(LatGPM)))

## Assign Date 
GPM.times <- paste0(substr(GPM.files,1,8),substr(GPM.files,11,14))
GPM.times [c(1,length(GPM.times))]
  
GPM.raster <- rasterize(GPM.coordinates2, MyGrid, values=as.array(t(GPM.list[["201610012100"]])))
  


##RAIN GUAGES-------------------------------------------------------------------

#L:È un altro campo georiferito di valori con un altro modo ancora di essere 
 #salvato. è un file binario
### Search DAT files 
RG.files <- dir(RG.folder, pattern=".dat")
RG.files <- RG.files[substr(RG.files,5,8)=="1001"]
length(RG.files)

### Assign Date 
RG.times <- substr(RG.files,1,12)
RG.times[c(1,length(RG.times)-1)]



### Import Pixel Coordinates (unstructured grid)
RG.coordinates <- read.table(file(paste0(datafolder, "PLV_coordinates.dat")),
                             sep = "", na.string = NA, as.is = TRUE, header = FALSE,
                             col.names = c("x","y","lat","lon"))

# trasformazioni necessarie ma ignote (non sappiamo l'origine dei dati)
# probabilmente sono stati salvati con una diversa convenzione di lettura di matrice
RG.coordinates$x <- rev(RG.coordinates$x)
RG.coordinates <- RG.coordinates %>% arrange(x, y)

### Import data (Precipitation Rate [mm/h])
nROW <- 444
nCOL <- 912
RG.df <- data.frame(matrix(ncol = 0, nrow = nROW*nCOL))
for(i in 1:(length(RG.files)) ) {
  tt <- RG.times[i]
  RG.df <- cbind(RG.df, readBin(paste0(RG.folder,RG.files[i]), "numeric", nROW*nCOL, size=4))  }

colnames(RG.df) <- RG.times

latlon <- tibble(id = 1:nrow(RG.coordinates), RG.coordinates[c("lon", "lat")])

latlon <- latlon %>% filter(
  lat >= 35,
  lat <= 50,
  lon >= 0,
  lon <= 25
)

RG.raster <- rasterize(cbind(latlon$lon, latlon$lat), MyGrid, values=RG.df[latlon$id,22])
RG.raster <- clamp(RG.raster, 0, values=F)
names(RG.raster) <- RG.times[22]
RG.raster
  #TABELLA

plot(RG.raster)
plot(ItalyBG, add =T) #, col = "grey85")

  #GRAFICO

#Dominio spaziale omogeneo
#L:Uso "mask" per ritagliare i prodotti che non mi interessano. Segue i bordi del
 #poligono
plot(mask(ERA5.MyGrid[["201610012100"]], ItalyBG))

ERA5.final <- mask(ERA5.MyGrid[["201610012100"]], ItalyBG)
RG.final <- mask(RG.raster, ItalyBG)

par(mfrow=c(1,2)) #L: Per fare combinazioni di plot
plot(ERA5.final)
plot(RG.final)
  









#ANDALISI PDF-------------------------------------------------------------------
#L:Devo creare prima ERA5.final
ERA5.data <- values(ERA5.final)
hist(ERA5.data)
  #ISTOGRAMMA
RG.data <- values(RG.final)
hist(RG.data)
#ISTOGRAMMA
tibble(
  Raingauges = RG.data,
  ERA5 = ERA5.data
) %>% 
  pivot_longer(
    cols = c("Raingauges", "ERA5"),
    names_to = "Product",
    values_to = "Precipitation"
  )%>%
  ggplot(aes(Precipitation, color = Product)) + 
  geom_freqpoly() + 
  scale_y_log10()

#ISTOGRAMMA




