# ---------------------------------------- #
# Create an animation of the rainfall maps #
# ---------------------------------------- #


# Clear the workspace
rm(list = ls())

# A new method to load packages:
# If they are not available in the system, they are automatically installed.

packages <- c(
  "tidyverse",
  "lubridate",
  "zoo",
  "terra",
  "rgdal",
  "ncdf4",
  "ggsci",
  "animation"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Define base folder paths
if(Sys.info()["sysname"]=="Linux"){
  rootfolder <- "/home/giacom0rovers1/tutorLabFisAtm22/Esp_3/"
  datafolder <- "~/Insync/giacomo.roversi2@studio.unibo.it/OneDrive Biz - Shared/Laboratorio di Fisica dell'Atmosfera/Studenti_Esp3/"
}else{
  rootfolder <- "C:/projects/tutorLabFisAtm22/Esp_3/"
  datafolder <- "C:/Users/grove/OneDrive - Alma Mater Studiorum Università di Bologna/Laboratorio di Fisica dell'Atmosfera/Studenti_Esp3/"
}

setwd(rootfolder)

resfolder <- "risultati/"
outfolder <- "extra/animazione/"

# Create the project directories, if missing
if (!dir.exists(outfolder)) {
  dir.create(outfolderfolder)
}

# defining mask
GlobalBG <- vect(paste0(datafolder, "World_Countries_(Generalized)"))
ItalyBG <- GlobalBG[112]

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

# read data
ERA5 <- rast(paste0(resfolder, "ERA5_final.tif"))
RG   <- rast(paste0(resfolder, "RG_final.tif"))
GPM  <- rast(paste0(resfolder, "GPM_final.tif"))
HSAF <- rast(paste0(resfolder, "HSAF_final.tif"))

# take only times with all data (HSAF has some holes)
names(RG)[which(!names(RG) %in% names(HSAF))]
new.times <- names(RG)[which(names(RG) %in% names(HSAF))]

RG   <- RG[[new.times]]
ERA5 <- ERA5[[new.times]]
GPM  <- GPM[[new.times]]


# plot function
plotMultiMap <- function(HSAF, GPM, ERA5, RG, layerNames){
  
  nCols <- 2
  nRows <- 2
  par(mfrow=c(nRows, nCols))
  
  # list with details for legend
  plg = list(
    title = expression(bold("Average rain\n rate (mm/h)")),
    title.cex = 0.7,
    cex = 0.7#,
    #shrink=0
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
        xlim = MyGrid_extLON,
        ylim = MyGrid_extLAT,
        range=c(0, 17),
        plg  = plg,
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

setwd(outfolder)
# make animation 
saveHTML({
  for(time in new.times[1:238]){
    plotMultiMap(HSAF, GPM, ERA5, RG, time)
  }
})