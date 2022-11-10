# --------------------------------------------------------------------------- #
#              Laboratorio di Fisica dell'Atmosfera - FST UniBo               #
#                     Data processing per l'esperienza 3                      #
#                 Giacomo Roversi (giacomo.roversi3@unibo.it)                 #
#                             10 novembre 2022                                #
# --------------------------------------------------------------------------- #

# Clear the workspace 
rm(list=ls())

# Timer
StartTime <- proc.time()

# ======================= CONFIG (parte da modificare) ========================

# Define base folder paths
rootfolder <- "/home/giacom0rovers1/tutorLabFisAtm22/Esp_3/"
datafolder <- "~/Insync/giacomo.roversi2@studio.unibo.it/OneDrive Biz - Shared/Laboratorio di Fisica dell'Atmosfera/Studenti_Esp3/"

# Define time interval (between 1 October and 30 November 2016)
sel_START <- "201610010000"
sel_END   <- "201610052300"

# Select properties of the common grid ("MyGrid"):
MyGrid_res = 0.1
MyGrid_extLON = c(0, 25)
MyGrid_extLAT = c(35, 50)

# Set the unic conversion coefficients
CorrectionFactor_HSAF = 3600 # kg/m2*s --> mm/h
CorrectionFactor_ERA5 = 1000 # m/h --> mm/h

# Select needed packages (default should work)
packages <- c("tidyverse",
              "lubridate",
              "zoo",
              "terra",
              "rgdal",
              "progress")

# ============================== FUNCTIONS ====================================

hourlyMean <- function(IN.raster, sel.times){
  #' Aggregates time layers of a SpatRaster object through an hourly mean.
  #' @param IN.raster SpatRaster object with multiple layers. Layer names should be in the "%Y%m%d%H%M" format.
  #' @param sel.times Array of strings of the desired time intervals (hours). This should also be in the "%Y%m%d%H%M" format.
  #' @return Output SpatRaster has the same properties of the input raster and layers corresponding to sel.times
  #' @author Giacomo Roversi
  
  # Read raster time intervals
  IN.times <- names(IN.raster)
  
  # Create empty objects
  OUT.raster <- rast(IN.raster)
  OUT.times  <- character() 
  
  # Progress bar
  bar <- progress_bar$new(format = "  :what [:bar] :current/:total (:percent) eta: :eta",
                          total = length(sel.times),
                          clear = FALSE)
  
  # Loop through requested time intervals
  for(tim in sel.times){
    
    # Create time strings for every minute in the selected hour
    interval <- format(seq(ymd_hm(tim) - hours(1) + minutes(1), 
                           ymd_hm(tim), 
                           by="1 min"), "%Y%m%d%H%M")
    
    # Select data rows within the selected hour
    idx <- which(IN.times %in% interval)
    
    if(length(idx) > 0){
      # Add the mean of the selected data to the output raster (new layer)
      OUT.raster <- c(OUT.raster, mean(IN.raster[[idx]]))
      OUT.times  <- c(OUT.times, tim)
    }
    
    bar$tick(tokens = list(what = "hourly mean "))
  }
  
  # Set layer names
  names(OUT.raster) <- OUT.times
  
  # Output
  return(OUT.raster) 
}


# ================================= INIT ======================================

# Install/load packages
for(pkg in packages){
  if(!require(pkg, character.only = T)){install.packages(pkg)} 
  library(pkg, character.only = T)
}

# Set the working directory
setwd(rootfolder)

# Data sources subdirectories (relative path)
HSAF.folder <- paste0(datafolder, "SAT/")
RG.folder   <- paste0(datafolder, "PLV/")
GPM.folder   <- paste0(datafolder, "GPM/")
ERA5.folder   <- paste0(datafolder, "ERA/")

# Project folders
figfolder  <- "figure/"
resfolder  <- "risultati/"

# Create the project directories, if missing
if(!dir.exists(figfolder)){  dir.create(figfolder) }
if(!dir.exists(resfolder)){  dir.create(resfolder) }


# Array of selected time intervals (hours)
sel.times <- format(seq(ymd_hm(sel_START)+hours(1), ymd_hm(sel_END), by="1 hour"), "%Y%m%d%H%M")

# Import the borders shapefiles of the countries of the world
# https://hub.arcgis.com/datasets/2b93b06dc0dc4e809d3c8db5cb96ba69_0/
GlobalBG <- vect(paste0(datafolder,"World_Countries_(Generalized)"))
ItalyBG  <- GlobalBG[112]

# Create the common grid
MyGrid <- rast(ncol=diff(MyGrid_extLON)/MyGrid_res, 
               xmin=MyGrid_extLON[1], xmax=MyGrid_extLON[2], 
               nrow=diff(MyGrid_extLAT)/MyGrid_res, 
               ymin=MyGrid_extLAT[1], ymax=MyGrid_extLAT[2])

# Format for progress bar inside rasterize for loops
barf <- "  :what [:bar] :current/:total (:percent) eta: :eta"

# ================================== RUN ======================================


# ERA5 ------------------------------------------------------------------------
cat(sprintf("ERA5\n"))

## Read file names ------------------------------------------------------------
ERA5.files <- dir(ERA5.folder, pattern=".grib")

## Import data ----------------------------------------------------------------
ERA5.raster <- rast(paste0(ERA5.folder, ERA5.files)) * CorrectionFactor_ERA5

## Assign Date ----------------------------------------------------------------
ERA5.times_START <- substr(ERA5.files,6,15)
ERA5.times_END   <- substr(ERA5.files,17,26)

ERA5.times <- seq(ymd_h(ERA5.times_START), 
                  ymd_h(ERA5.times_END), 
                  by="1 hour") %>% format("%Y%m%d%H%M")
names(ERA5.raster) <- ERA5.times

## Project over MyGrid  -------------------------------------------------------
ERA5.MyGrid <- terra::project(ERA5.raster, MyGrid)

## Filter out negative values -------------------------------------------------
ERA5.MyGrid <- clamp(ERA5.MyGrid, 0, values=F)

## Select only Italy and desired time intervals  ------------------------------
ERA5.final <- mask(ERA5.MyGrid[[sel.times]], ItalyBG)

## Save GeoTiff file ----------------------------------------------------------
writeRaster(ERA5.final, paste0(resfolder, "ERA5_final.tif"), overwrite=TRUE )

## Clear memory  ---------------------------------------------------------------
rm(ERA5.raster, ERA5.MyGrid, ERA5.final)




# H-SAF -----------------------------------------------------------------------
cat(sprintf("H-SAF\n"))

## Search GRIB files ----------------------------------------------------------
HSAF.files <- dir(HSAF.folder, pattern=".grb")

## Assign Date ----------------------------------------------------------------
HSAF.times <- paste0(substr(HSAF.files,5,12),substr(HSAF.files,14,17))

## Select only files inside the selected time interval -------------------------
sel_idx <- which(ymd_hm(HSAF.times) >= ymd_hm(sel_START) & 
                   ymd_hm(HSAF.times) <= ymd_hm(sel_END))
HSAF.files <- HSAF.files[sel_idx]
HSAF.times <- HSAF.times[sel_idx]

## Import Coordinates (unstructured grid) -------------------------------------
HSAF.coordinates <- read.table(file(paste0(datafolder,"SAT_coordinates.dat")),
                               sep = " ", 
                               na.string = NA, 
                               as.is = TRUE, 
                               header = FALSE)

HSAF.df <- HSAF.coordinates[1710000:1,]
bar <- progress_bar$new(format = barf, 
                        total = length(HSAF.times),
                        clear = FALSE)

for(i in 1:length(HSAF.times)){
  HSAF.df <- cbind(HSAF.df, 
                   matrix(rgdal::readGDAL(paste0(HSAF.folder,HSAF.files[i]), 
                                          silent = TRUE)[[1]], 1710000, 1))
  bar$tick(tokens = list(what = "reading data"))
} 
colnames(HSAF.df) <- c("Lon", "Lat", HSAF.times)

## Restrict domain to Mediterranean area (bulk) -------------------------------
HSAF.df_ITA <- HSAF.df %>% filter(
  Lon >= 0,
  Lon <= 30,
  Lat >= 30,
  Lat <= 60
) 

## Create SpatRaster object from MyGrid ---------------------------------------
HSAF.raster <- rast(MyGrid)

## Rasterize data over SpatRaster object (layr by layer) ----------------------
bar <- progress_bar$new(format = barf, 
                        total = length(HSAF.times),
                        clear = FALSE)

for(i in 1:length(HSAF.times)){
  HSAF.raster <- c(HSAF.raster, 
                   rasterize(cbind(HSAF.df_ITA$Lon, HSAF.df_ITA$Lat), 
                             MyGrid, 
                             values=as.array(HSAF.df_ITA[[HSAF.times[i]]])
                   ) * CorrectionFactor_HSAF)
  bar$tick(tokens = list(what = "rasterize   "))
}

## Set layer names ------------------------------------------------------------
names(HSAF.raster) <- HSAF.times

## Filter out negative values -------------------------------------------------
HSAF.raster <- clamp(HSAF.raster, 0, values=F)

## Aggregate to hourly averages -----------------------------------------------
HSAF.hourly <- hourlyMean(HSAF.raster, sel.times)

## Select only pixels inside Italy --------------------------------------------
HSAF.final  <- mask(HSAF.hourly, ItalyBG)

## Save GeoTiff file ----------------------------------------------------------
writeRaster(HSAF.final, paste0(resfolder, "HSAF_final.tif"), overwrite=TRUE )

## Clear memory  --------------------------------------------------------------
rm(HSAF.df, HSAF.df_ITA, HSAF.raster, HSAF.hourly, HSAF.final, HSAF.coordinates)




# GPM-IMERG -------------------------------------------------------------------
cat(sprintf("GPM-IMERG\n"))

## Search IMR files (ASCII) ---------------------------------------------------
GPM.files <- dir(GPM.folder, pattern=".imr")

## Assign Date ----------------------------------------------------------------
GPM.times <- paste0(substr(GPM.files,1,8),substr(GPM.files,11,14))

## Select only files inside the selected time interval -------------------------
sel_idx <- which(ymd_hm(GPM.times) >= ymd_hm(sel_START) & 
                   ymd_hm(GPM.times) <= ymd_hm(sel_END))
GPM.files <- GPM.files[sel_idx]
GPM.times <- GPM.times[sel_idx]

## Import data (Precipitation Rate [mm/h]) ------------------------------------
GPM.list <- list()
bar <- progress_bar$new(format = barf, 
                        total = length(GPM.files),
                        clear = FALSE)

for(i in 1:(length(GPM.files)) ) {
  GPM.list[[i]] <- read.table(file(paste0(GPM.folder,GPM.files[i])), 
                              skip = 181, 
                              nrows = 180,
                              sep = ",", 
                              na.string = NA, 
                              as.is = TRUE, 
                              header = FALSE)
  
  GPM.list[[i]]  <- as.matrix(GPM.list[[i]][-1])
  bar$tick(tokens = list(what = "reading data"))
}
names(GPM.list) <- GPM.times

## Import Coordinates (structured grid) ---------------------------------------
LatGPM <- read.table(file(paste0(GPM.folder,GPM.files[1])), 
                     skip = 361, 
                     nrows = 1,
                     sep = ",", 
                     na.string = NA, 
                     as.is = TRUE, 
                     header = FALSE)
temp   <- LatGPM[[2]]
for(i in 3:length(LatGPM)) temp[i-1] <- LatGPM[[i]]
LatGPM <- temp

LonGPM <- read.table(file(paste0(GPM.folder,GPM.files[1])), 
                     skip = 362, 
                     nrows = 1,
                     sep = ",", 
                     na.string = NA, 
                     as.is = TRUE, 
                     header = FALSE)
temp   <- LonGPM[[2]]
for(i in 3:length(LonGPM)) temp[i-1] <- LonGPM[[i]]
LonGPM <- temp

GPM.coordinates <- 0
for(i in LonGPM) for(j in LatGPM) GPM.coordinates <- c(GPM.coordinates,i,j)
GPM.coordinates <- t(matrix(GPM.coordinates[-1], 2, length(LonGPM)*length(LatGPM)))

## Create SpatRaster object from MyGrid ---------------------------------------
GPM.raster <- rast(MyGrid)

## Rasterize data over SpatRaster object (layr by layer) ----------------------
bar <- progress_bar$new(format = barf,
                        total = length(GPM.times),
                        clear = FALSE)

for(i in 1:length(GPM.times)){
  GPM.raster <- c(GPM.raster, 
                  rasterize(GPM.coordinates, 
                            MyGrid, 
                            values=as.array(t(GPM.list[[GPM.times[i]]]))))
  bar$tick(tokens = list(what = "rasterize   "))
}
## Set layer names ------------------------------------------------------------
names(GPM.raster) <- GPM.times

## Filter out negative values -------------------------------------------------
GPM.raster <- clamp(GPM.raster, 0, values=F)

## Aggregate to hourly averages -----------------------------------------------
GPM.hourly <- hourlyMean(GPM.raster, sel.times)

## Select only pixels inside Italy --------------------------------------------
GPM.final  <- mask(GPM.hourly, ItalyBG)

## Save GeoTiff file ----------------------------------------------------------
writeRaster(GPM.final, paste0(resfolder, "GPM_final.tif"), overwrite=TRUE )

## Clear memory  --------------------------------------------------------------
rm(GPM.list, GPM.raster, GPM.hourly, GPM.final, GPM.coordinates, LatGPM, LonGPM)



# RAIN GAUGES -----------------------------------------------------------------
cat(sprintf("RainGauges\n"))

## Search DAT files -----------------------------------------------------------
RG.files <- dir(RG.folder, pattern=".dat")

## Assign Date ----------------------------------------------------------------
RG.times <- substr(RG.files,1,12)

## Select only files inside the selected time interval -------------------------
sel_idx <- which(ymd_hm(RG.times) >= ymd_hm(sel_START) & 
                   ymd_hm(RG.times) <= ymd_hm(sel_END))
RG.files <- RG.files[sel_idx]
RG.times <- RG.times[sel_idx]

## Import Pixel Coordinates (unstructured grid) -------------------------------
RG.coordinates <- read.table(file(paste0(datafolder, "PLV_coordinates.dat")),
                             sep = "", na.string = NA, as.is = TRUE, header = FALSE,
                             col.names = c("x","y","lat","lon"))

# Transformations for consistent data indexing
RG.coordinates$x <- rev(RG.coordinates$x)
RG.coordinates <- RG.coordinates %>% arrange(x, y)

## Import data (Precipitation Rate [mm/h]) ------------------------------------
nROW <- 444
nCOL <- 912
RG.df <- data.frame(matrix(ncol = 0, nrow = nROW*nCOL))

bar <- progress_bar$new(format = barf, 
                        total = length(RG.files),
                        clear = FALSE)

for(i in 1:(length(RG.files)) ) {
  tt <- RG.times[i]
  RG.df <- cbind(RG.df, readBin(paste0(RG.folder,RG.files[i]), 
                                "numeric", 
                                nROW*nCOL, 
                                size=4))
  bar$tick(tokens = list(what = "reading data"))
}
colnames(RG.df) <- RG.times

# Create support data frame with coordinates and id (used for filtering)
latlon <- tibble(id = 1:nrow(RG.coordinates), RG.coordinates[c("lon", "lat")])

## Restrict domain to Mediterranean area (bulk) -------------------------------
latlon <- latlon %>% filter(
  lat >= 35,
  lat <= 50,
  lon >= 0,
  lon <= 25
)

## Create SpatRaster object from MyGrid ---------------------------------------
RG.raster <- rast(MyGrid)

## Rasterize data over SpatRaster object (layr by layer) ----------------------
bar <- progress_bar$new(format = barf,
                        total = length(RG.times),
                        clear = FALSE)

for(i in 1:length(RG.times)){
  RG.raster <- c(RG.raster, 
                 rasterize(cbind(latlon$lon, latlon$lat), 
                           MyGrid, 
                           values=RG.df[latlon$id,22]))
  bar$tick(tokens = list(what = "rasterize   "))
}

## Set layer names ------------------------------------------------------------
names(RG.raster) <- RG.times

## Filter out negative values -------------------------------------------------
RG.raster <- clamp(RG.raster, 0, values=F)

## Select only Italy and desired time intervals  ------------------------------
RG.final   <- mask(RG.raster[[sel.times]], ItalyBG)

## Save GeoTiff file ----------------------------------------------------------
writeRaster(RG.final, paste0(resfolder, "RG_final.tif"), overwrite=TRUE )

## Clear memory  --------------------------------------------------------------
rm(RG.df, RG.raster, RG.final, RG.coordinates, latlon)


# All done --------------------------------------------------------------------
cat(sprintf("All done. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))
