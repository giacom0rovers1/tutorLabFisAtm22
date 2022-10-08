# ===========================================================
# = LABORATORY OF DATA PROCESSING AND EVALUATION            =
# = Precipitation from several data sources     v2020.11.09 =
# ===========================================================


# Needs rGDAL Package (Geospatial Data Abstraction Library)
# install.packages(rgdal)     # Just at the first time
# install.packages(gdalUtilities) # Just at the first time

rm(list=ls(all=TRUE))
gc(reset=TRUE)
library(rgdal)      # Allows call a command directly with no need to indicate a library
library(gdalUtilities)  # (library::command, only commands out of the R base)

# Set Working Directory --------------------------------------------------------
setwd("C:/Users/qucci/Desktop/RTD_2017/Corso Lab/Corso 2021_22")

# Background Maps --------------------------------------------------------------
ContoursGLOBAL <- dir("C:/Users/qucci/Desktop/RTD_2017/Corso Lab/Corso 2021_22/Contours Global")
ContoursIT     <- dir("C:/Users/qucci/Desktop/RTD_2017/Corso Lab/Corso 2021_22/Contours IT")

GlobalBG <- function(ColFILL,ColEDGE)
              for (i in 1:length(ContoursGLOBAL)) {
                temp <- read.table(paste0("./Contours Global/",ContoursGLOBAL[i]))
                polygon(x = temp[,1], y = temp[,2], lwd = 1,
                col = ColFILL, border = ColEDGE)  }

ItalyBG  <- function(ColFILL,ColEDGE)
              for (i in 1:length(ContoursIT)) {
                temp <- read.table(paste0("./Contours IT/",ContoursIT[i]))
                polygon(x = temp[,1], y = temp[,2], lwd = 1,
                col = ColFILL, border = ColEDGE)  }

# Legends to Precipitation -----------------------------------------------------
PrecLEG <- function(precUNIT)
           legend("bottomright",
             c("<  2", "<  4", "<  8", "< 16", "> 16"),
             col = c("blue","green","yellow","orange","red"),
             bg = rgb(.95,.95,.95), pch = 15, border = NA, cex = 0.75,
             title = paste0("[",precUNIT,"]")  )

asDATE <- function(tt) paste0(substr(tt,7,8),"/",substr(tt,5,6),"/",
                         substr(tt,1,4),"  ",substr(tt,9,10),":",
                         substr(tt,11,12),"Z")



# A) PRECIPITATION FROM SATELLITE DATA (H-SAF) =================================

# Search GRIB files ------------------------------------------------------------
fSAT <- dir("SAT", pattern=".grb")
fSAT <- fSAT[substr(fSAT,9,12)=="1001"]
length(fSAT)

# Import data (Precipitation Rate [kg/m?s]) ------------------------------------
pSAT <- list()
for(i in 1:length(fSAT)) 
  pSAT[[i]] <- readGDAL(paste0("SAT/",fSAT[i]), silent = TRUE)

length(pSAT)
dim(pSAT[[1]])
image(pSAT[[1]])

# Import Coordinates (unstructured grid) ---------------------------------------
cSAT <- read.table(file("SAT_coordinates.dat"),
                   sep = " ", na.string = NA, as.is = TRUE, header = FALSE)
dim(cSAT)

# Assign Date ------------------------------------------------------------------
tSAT <- paste0(substr(paste0("SAT/",fSAT),9,16),substr(paste0("SAT/",fSAT),18,21))
tSAT [c(1,length(tSAT))]

# Merge Data -------------------------------------------------------------------
temp <- cSAT[1710000:1,]
for(i in 1:length(tSAT)) temp <- cbind(temp,matrix(pSAT[[i]][[1]],dim(pSAT[[1]])[1],1))
pSAT <- temp
colnames(pSAT) <- c("Lon","Lat",tSAT)

pSAT[1:8,1:6]
dim(pSAT)

# Domain Reduction to Italy ----------------------------------------------------
pSAT <- pSAT[(pSAT[,1]>= 0)&(pSAT[,1]<=30)&(pSAT[,2]>=30)&(pSAT[,2]<=60),]
rownames(pSAT) <- 1:dim(pSAT)[1]
pSAT[1:8,1:5]
dim(pSAT)

rm(list=c("cSAT","fSAT","tSAT","temp"))   # Clear useless memory
gc(reset=TRUE)

# Precipitation Rate Plot ------------------------------------------------------
PlotSAT <- function(tt,minLON,maxLON,minLAT,maxLAT,BG) {tt <- tt+2
           plot(c(0,0), type = "n", asp = 1.2, cex.main = 0.75,
                xlim = c(minLON,maxLON), ylim = c(minLAT,maxLAT),
                xlab = "Longitude (°E)", ylab = "Latitude (°N)",
                main = paste0("PRECIPITATION FROM H-SAF DATA (SATELLITE)\n",
                  substr(colnames(pSAT)[tt],7,8),"/",substr(colnames(pSAT)[tt],5,6),"/",
                  substr(colnames(pSAT)[tt],1,4),"  ",substr(colnames(pSAT)[tt],9,10),":",
                  substr(colnames(pSAT)[tt],11,12),"Z"),
                panel.first = {
                  BG(rgb(.95,.95,.95),"gray")
                  grid()
                  for(i in 1:5)
                    points(pSAT[(pSAT[,tt]> 0.0001*c(.1,2,4,8,16,1000)[i])&
                                (pSAT[,tt]<=0.0001*c(.1,2,4,8,16,1000)[i+1]),c(1,2)],
                           pch = 19, cex = 0.05, col = c("blue","green","yellow","orange","red")[i])
                  PrecLEG("x10-4 kg/m?s") } )}

#expression(paste("x10-4 kg/m" ^"2", "s", sep = " "))

PlotSAT(5,0,30,30,60,GlobalBG)

# B) PRECIPITATION FROM PLUVIOMETRIC DATA (KRIGING INTERPOLATION) ==============

# Search DAT files -------------------------------------------------------------
fPLV <- dir("PLV", pattern=".dat")
fPLV <- fPLV[substr(fPLV,5,8)=="1001"]
length(fPLV)

# Import data (Precipitation Rate [mm/h]) --------------------------------------
pPLV <- list()
nROW <- 444
nCOL <- 912
for(i in 1:(length(fPLV)) ) {
  pPLV[[i]] <- readBin(paste0("PLV/",fPLV[i]), "numeric", nROW*nCOL, size = 4)  }

# Import Pixel Coordinates (unstructured grid) ---------------------------------
cPLV <- read.table(file("PLV_coordinates.dat"),
                   sep = "", na.string = NA, as.is = TRUE, header = FALSE,
                   col.names = c("x","y","lat","lon")
                   )
cPLV$x <- rev(cPLV$x)

# Assign Date ------------------------------------------------------------------
tPLV <- substr(paste0("PLV/",fPLV),5,16)
tPLV[c(1,length(tPLV)-1)]

# Pixel Borders Plot -----------------------------------------------------------
PixCOR  <- function(pX,pY)
           dplyr::filter(cPLV,(cPLV$x >= pX-1)&(cPLV$x <= pX)&
                              (cPLV$y >= pY-1)&(cPLV$y <= pY))[c(1,2,4,3),c(4,3)]
PixBOR  <- function(pX,pY) {
           polygon(PixCOR(pX,pY), lwd=2, col = NA, border = "lightblue")
           points(PixCOR(pX,pY), pch = 20, col = "blue")
           points(mean(PixCOR(pX,pY)[,1]),mean(PixCOR(pX,pY)[,2]),
           pch = 13, col = "red") }

plot(PixCOR(1,1),type = "n", panel.first = c( grid(),PixBOR(1,1) )  )

# Precipitation as Pixel Plot --------------------------------------------------
PixCONT <- function(tt,pID) if(pPLV[[tt]][pID] > 0){
             pX  <- floor(pID/nCOL)+1
             pY  <- (pID %% nCOL)
             pXY <- dplyr::filter(cPLV,(cPLV$x >= pX-1)&(cPLV$x <= pX)&
                                       (cPLV$y >= pY-1)&(cPLV$y <= pY))
             pXY <- pXY[c(1,2,4,3),c(4,3)]
           # points(mean(pXY[,1]),mean(pXY[,2]), pch = 19, col = "red", cex = 0.3)
             polygon(pXY, lwd=2, border = NA,
                     col = ifelse(pPLV[[tt]][pID] <= 2, rgb(0,0,1,.5),
                           ifelse(pPLV[[tt]][pID] <= 4, rgb(0,1,0,.5),
                           ifelse(pPLV[[tt]][pID] <= 8, rgb(1,1,0,.5),
                           ifelse(pPLV[[tt]][pID] <= 16,rgb(1,.65,0,.5),
                                                        rgb(1,0,0,.5) ))))  )}

PlotPLV <- function (tt,minLON,maxLON,minLAT,maxLAT,BG)
           plot(c(0,0), type = "n", asp = 1.2, cex.main = 0.75,
                xlim = c(minLON,maxLON), ylim = c(minLAT,maxLAT),
                xlab = "Longitude (°E)", ylab = "Latitude (°N)",
                main = paste0("PRECIPITATION FROM PLUVIOMETRIC DATA\n",
                         substr(tPLV[tt],7,8),"/",substr(tPLV[tt],5,6),"/",
                         substr(tPLV[tt],1,4),"  ",substr(tPLV[tt],9,10),":",
                         substr(tPLV[tt],11,12),"Z"),
                panel.first = {BG(rgb(.95,.95,.95),"gray")
                               grid()
                               for(i in 1:length(pPLV[[1]])) PixCONT(tt,i)
                               PrecLEG("mm/h")
                               })
PlotPLV(22,0,30,30,60,GlobalBG)



# C) PRECIPITATION FROM GPM DATA (SATELLITE) ===================================

# Search IMR files (ASCII) -----------------------------------------------------
fGPM <- dir("GPM", pattern=".imr")
fGPM <- fGPM[substr(fGPM,5,8)=="1001"]
length(fGPM)

# Import data (Precipitation Rate [mm/h]) --------------------------------------
pGPM <- list()
for(i in 1:(length(fGPM)) ) {
  pGPM[[i]] <- read.table(file(paste0("GPM/",fGPM[i])), skip = 181, nrows = 180,
                    sep = ",", na.string = NA, as.is = TRUE, header = FALSE)
  pGPM[[i]]  <- as.matrix(pGPM[[i]][-1])
}

length(pGPM)
dim(pGPM[[1]])

# Import Coordinates (structured grid) -----------------------------------------
LatGPM <- read.table(file(paste0("GPM/",fGPM[1])), skip = 361, nrows = 1,
                     sep = ",", na.string = NA, as.is = TRUE, header = FALSE)
temp   <- LatGPM[[2]]
for(i in 3:length(LatGPM)) temp[i-1] <- LatGPM[[i]]
LatGPM <- temp

LonGPM <- read.table(file(paste0("GPM/",fGPM[1])), skip = 362, nrows = 1,
                     sep = ",", na.string = NA, as.is = TRUE, header = FALSE)
temp   <- LonGPM[[2]]
for(i in 3:length(LonGPM)) temp[i-1] <- LonGPM[[i]]
LonGPM <- temp

c(length(LonGPM),length(LatGPM))

cGPM <- 0
for(i in LonGPM) for(j in LatGPM) cGPM <- c(cGPM,i,j)
cGPM <- t(matrix(cGPM[-1], 2, length(LonGPM)*length(LatGPM)))

# Assign Date ------------------------------------------------------------------
tGPM <- paste0(substr(paste0("GPM/",fGPM),5,12),substr(paste0("GPM/",fGPM),15,19))
tGPM [c(1,length(tGPM))]


# Precipitation Rate Plot ------------------------------------------------------

PlotGPM <- function(tt,minLON,maxLON,minLAT,maxLAT,BG) {
           plot(c(0,0), asp = 1.2, cex.main = 0.75,
             xlim = c(minLON,maxLON), ylim = c(minLAT,maxLAT),
             xlab = "Longitude (°E)", ylab = "Latitude (°N)",
             main = paste0("PRECIPITATION FROM GPM DATA (SATELLITE)\n",
                           asDATE(tGPM[tt])),
             panel.first = BG(rgb(.95,.95,.95),"gray"))
           .filled.contour(LonGPM, LatGPM, pGPM[[tt]],
             levels = c(1,2,4,8,16,1000),
             col = c(rgb(0,0,1,.75),rgb(0,1,0,.75),rgb(1,1,0,.75),
                     rgb(1,.65,0,.75),rgb(1,0,0,.75) )  )
           grid()
           PrecLEG("mm/h")    }
PlotGPM(45,0,30,30,60,GlobalBG)


# D) PRECIPITATION FROM ERA5 DATA ==============================================

# Search GRIB files ------------------------------------------------------------
fERA5 <- dir(path = "C:/Users/qucci/Desktop/RTD_2017/Corso Lab/Corso 2021_22", pattern=".grib")
length(fERA5)

# Import data ------------------------------------------------------------------
pERA5 <- readGDAL(fERA5, silent = TRUE)
dim(pERA5)

# Convert units from 'm' to 'mm' -----------------------------------------------
for (i in 1:dim(pERA5)[2]) { pERA5[[i]] <- 1000*pERA5[[i]] }

# Assign Date ------------------------------------------------------------------
tERA5        <- as.numeric(substr(fERA5,6,15))
for(h in 1:dim(pERA5)[2]-1)
  tERA5[h+1] <- tERA5[1] + h%%24 + 100*floor(h/25) + ifelse(h%%24==0,100,0)
tERA5        <- (tERA5-100)*100
names(pERA5) <- tERA5

# Precipitation as Contour Plot ------------------------------------------------
PlotERA5 <- function(tt,minLON,maxLON,minLAT,maxLAT,BG) {
            xLon   <- unique(coordinates(pERA5)[,1])
            yLat   <- unique(rev(coordinates(pERA5)[,2]))
            ppERA5 <- matrix(pERA5[[tt]],length(xLon),length(yLat))
            temp   <- matrix(pERA5[[tt]],length(xLon),length(yLat))
            for(i in 1:length(xLon)) for(j in 1:length(yLat))
            ppERA5[i,length(yLat)-j+1] <- temp[i,j]
            plot(c(0,0), type = "n", asp = 1.2, cex.main = 0.75,
              xlim = c(minLON,maxLON), ylim = c(minLAT,maxLAT),
              xlab = "Longitude (°E)", ylab = "Latitude (°N)",
              main = paste0("PRECIPITATION FROM ERA5 DATA\n",
                            asDATE(tERA5[tt])),
              panel.first = BG(rgb(.95,.95,.95),"gray")  )
            .filled.contour(xLon, yLat, ppERA5,
              levels = c(1,2,4,8,16,1000),
              c(rgb(0,0,1,.5),rgb(0,1,0,.5),rgb(1,1,0,.5),
                rgb(1,.65,0,.5),rgb(1,0,0,.5)))
            grid()
            PrecLEG("mm/h") }
PlotERA5(22,0,30,30,60,GlobalBG)


# DOMAIN COMPARISON ============================================================

# Domain Dimension -------------------------------------------------------------
plot(c(0,0), type = "n", asp = 1.2, cex.main = 0.75,
     xlim = c(-15,55), ylim = c(28,62),
     xlab = "Longitude (?E)", ylab = "Latitude (?N)",
     main = "PRECIPITATION DOMAIN COMPARISON\nand number of gridpoints/pixels",
     panel.first = {GlobalBG(rgb(.95,.95,.95),"gray")
                    grid()  })

# A) H-SAF DATA (SATELLITE) DOMAIN
     polygon(rbind(c(min(pSAT[,1]),min(pSAT[,2])),
       c(max(pSAT[,1]),min(pSAT[,2])),
       c(max(pSAT[,1]),max(pSAT[,2])),
       c(min(pSAT[,1]),max(pSAT[,2]))),
       lwd = 2, col = NA, border = "orange")
     text(27, 61.5, "H-SAF",      cex = 0.75, col = "orange")
     text(26, 28.5, dim(pSAT)[1], cex = 0.75, col = "orange")

# B) PLUVIOMETRIC DATA DOMAIN
     polygon(rbind(cPLV[(cPLV[,2]%%nCOL==0),c(4,3)],
       cPLV[(cPLV[,2]==(nCOL-1)),c(4,3)][nROW:1,]),
     lwd = 2, col = NA, border = "red")
     text(48,54,"Pluviometric",     cex = 0.75, col = "red")
     text(41,34, length(pPLV[[1]]), cex = 0.75, col = "red")

# C) GPM DATA (SATELLITE) DOMAIN
     polygon(rbind(c(min(LonGPM),min(LatGPM)),
       c(max(LonGPM),min(LatGPM)),
       c(max(LonGPM),max(LatGPM)),
       c(min(LonGPM),max(LatGPM))),
       lwd = 2, col = NA, border = "brown")
     text(20,50.2, "GPM",             cex = 0.75, col = "brown")
     text(19,  33, length(pGPM[[1]]), cex = 0.75, col = "brown")

# D) ERA5 DATA DOMAIN
     polygon(rbind(
       c(min(coordinates(pERA5)[,1]),min(coordinates(pERA5)[,2])),
       c(max(coordinates(pERA5)[,1]),min(coordinates(pERA5)[,2])),
       c(max(coordinates(pERA5)[,1]),max(coordinates(pERA5)[,2])),
       c(min(coordinates(pERA5)[,1]),max(coordinates(pERA5)[,2]))),
       lwd = 2, col = NA, border = "blue")
     text(15.8,47,"ERA5", cex = 0.75, col = "blue")
     text(9,37.5,length(pERA5[[1]]), cex = 0.75, col = "blue")

# Domain Gridpoints (La Spezia Case) -------------------------------------------
plot(c(0,0), type = "n", asp = 1.2, cex.main = 0.75,
     xlim = c(9.6,10.2), ylim = c(43.9,44.35),
     xlab = "Longitude (°E)", ylab = "Latitude (°N)",
     main = paste0("PRECIPITATION GRIDPOINTS COMPARISON\n",
             "La Spezia (Italy) 9.6-10.2°E and 43.9-44.3°N"),
     panel.first = {ItalyBG(rgb(.95,.95,.95),"gray")
                    grid()  })

    # A) H-SAF DATA (SATELLITE) GRIDPOINTS 
         points(pSAT[(pSAT[,1]>=9.6)&(pSAT[,1]<=10.20)&
                     (pSAT[,2]>=43.90)&(pSAT[,2]<=44.30),c(1,2)],
                      pch = 19, col = "orange", cex = 0.5)

    # B) PLUVIOMETRIC DATA GRIDPOINTS 
         points(cPLV[(cPLV$lon>=9.6)&(cPLV$lon<=10.20)&
                     (cPLV$lat>=43.90)&(cPLV$lat<=44.30),c(4,3)],
                     pch = 19, col = "red", cex = 0.5)

    # C) GPM DATA (SATELLITE) GRIDPOINTS 
         points(cGPM[(cGPM[,1]>=9.6)&(cGPM[,1]<=10.20)&
                     (cGPM[,2]>=43.90)&(cGPM[,2]<=44.30),],
                pch = 19, col = "brown", cex = 0.5)

    # D) ERA5 DATA GRIDPOINTS 
         points(coordinates(pERA5)[(coordinates(pERA5)[,1]>=9.6)&
                                   (coordinates(pERA5)[,1]<=10.20)&
                                   (coordinates(pERA5)[,2]>=43.90)&
                                   (coordinates(pERA5)[,2]<=44.30),],
                                   pch = 19, col = "blue", cex = 0.75)

    # E) NUMBER OF GRIDPOINTS
         text( 9.81, 44.33, paste0("H-SAF (",
              length(pSAT[(pSAT[,1]>=9.6)&(pSAT[,1]<=10.20)&
                          (pSAT[,2]>=43.90)&(pSAT[,2]<=44.30),1]),")"),
              cex = 0.75, col = "orange")
         text( 9.93, 44.33, paste0("PLV (",
              length(cPLV[(cPLV$lon>=9.6)&(cPLV$lon<=10.20)&
                          (cPLV$lat>=43.90)&(cPLV$lat<=44.30),1]),")"),
              cex = 0.75, col = "red")
         text(10.05, 44.33, paste0("GPM (",
              length(LonGPM[(LonGPM>=9.6)&(LonGPM<=10.20)])*
              length(LatGPM[(LatGPM>=43.90)&(LatGPM<=44.30)]),")"),
              cex = 0.75, col = "brown")
         text(10.17, 44.33, paste0("ERA5 (",
              length(coordinates(pERA5)[
                    (coordinates(pERA5)[,1]>=9.6)&
                    (coordinates(pERA5)[,1]<=10.20)&
                    (coordinates(pERA5)[,2]>=43.90)&
                    (coordinates(pERA5)[,2]<=44.30),1]),")"),
              cex = 0.75, col = "blue")



# COMMON GRID OF 0.1? USING NEAREST NEIGHBOUR (LA SPEZIA) ======================

# Build La Spezia's grid -------------------------------------------------------
MyGrid <- 0.1
for(i in seq(9.6,10.2, MyGrid[1]))
  for(j in seq(43.9,44.3, MyGrid[1]))
    MyGrid <- c(MyGrid,i,j)
MyGrid <- as.matrix(t(matrix(MyGrid[-1],2,length(MyGrid[-1])/2)))
rownames(MyGrid) <- 1:dim(MyGrid)[1]
colnames(MyGrid) <- c("lon","lat")

points(MyGrid, pch = 3)


# Define the Haversine function (distance between 2 points on a sphere) --------
dHAV <- function(lon1,lat1,lon2,lat2)
        6371000*sqrt((  ( ((lon2*pi/180)-(lon1*pi/180))*
                     cos( ((lat1*pi/180)+(lat2*pi/180) )/2)  )^2 +
                         ( (lat2*pi/180)-(lat1*pi/180) )^2))

# Define a function to find the nearest neighbour ------------------------------
NearGP <- function(iMyGrid,Coords) {
          nearGP  <- c(-180,-90)
          mindist <- 999999
          temp <- Coords[(Coords[,1]>=MyGrid[iMyGrid,1]-.5)&
                         (Coords[,1]<=MyGrid[iMyGrid,1]+.5)&
                         (Coords[,2]>=MyGrid[iMyGrid,2]-.5)&
                         (Coords[,2]<=MyGrid[iMyGrid,2]+.5),]
          for(i in 1:dim(temp)[1])
           if(dHAV(MyGrid[iMyGrid,1],MyGrid[iMyGrid,2],temp[i,1],temp[i,2]) < mindist){
             nearGP  <- temp[i,]
             mindist <- dHAV(MyGrid[iMyGrid,1],MyGrid[iMyGrid,2],temp[i,1],temp[i,2]) }
          temp <- matrix(c(nearGP,mindist),1,3)
          colnames(temp) <- c("   Lon(°E)","   Lat(°N)","   Dist(m)")
          temp}

points(MyGrid[13,1],MyGrid[13,2], pch = 19)

NearGP(13,pSAT[,c(1,2)])
lines(rbind(MyGrid[13,],NearGP(13,pSAT[,c(1,2)])[,c(1,2)]), col="orange")

NearGP(13,cPLV[,c(4,3)])
lines(rbind(MyGrid[13,],NearGP(13,cPLV[,c(4,3)])[,c(1,2)]), col="red")

NearGP(13,cGPM)
lines(rbind(MyGrid[13,],NearGP(13,cGPM)[,c(1,2)]), col="brown")

NearGP(13,coordinates(pERA5))
lines(rbind(MyGrid[13,],NearGP(13,coordinates(pERA5))[,c(1,2)]), col="blue")

points(MyGrid[13,1],MyGrid[13,2], pch = 19)

# Build MyGrid from several datasets -------------------------------------------

# A) H-SAF DATA (SATELLITE) GRIDPOINTS ---------- 
MyGridSAT <- 0
for(tt in 3:length(pSAT))
  for(i in 1:dim(MyGrid)[1])
    MyGridSAT <- c(MyGridSAT,pSAT[ (pSAT[,1]==NearGP(i,pSAT[,c(1,2)])[[1]])&
                                   (pSAT[,2]==NearGP(i,pSAT[,c(1,2)])[[2]]), tt])

MyGridSAT <- matrix(MyGridSAT[-1],dim(MyGrid)[1])
colnames(MyGridSAT) <- colnames(pSAT)[-c(1,2)]
rownames(MyGridSAT) <- 1:dim(MyGrid)[1]
dim(MyGridSAT)

# B) PLUVIOMETRIC DATA GRIDPOINTS ---------------
cPixPLV <- read.table("PLV_central_coordinates.dat", sep = "\t")
dim(cPixPLV)

MyGridPLV <- 0
for(tt in 1:length(pPLV))
  for(i in 1:dim(MyGrid)[1])
    MyGridPLV <- c(MyGridPLV,
                   pPLV[[tt]][
                     as.numeric(sort(rownames(cPixPLV[cPixPLV[,1]==
                     NearGP(i,cPixPLV[(is.na(cPixPLV[,1])==FALSE)&
                                      (is.na(cPixPLV[,2])==FALSE),] )[[1]],]
                     ))[1])   ])

MyGridPLV <- matrix(MyGridPLV[-1],dim(MyGrid)[1])
colnames(MyGridPLV) <- tPLV
rownames(MyGridPLV) <- 1:dim(MyGrid)[1]
dim(MyGridPLV)

# C) GPM DATA (SATELLITE) GRIDPOINTS ------------
LonGPM2 <- cbind(LonGPM,1:length(LonGPM))
LatGPM2 <- cbind(LatGPM,1:length(LatGPM))

MyGridGPM <- 0
for(tt in 1:length(pGPM))
  for(i in 1:dim(MyGrid)[1])
    MyGridGPM <- c(MyGridGPM, pGPM[[tt]][ 
                   LonGPM2[LonGPM2[,1]==NearGP(i,cGPM[,c(1,2)])[[1]],2],
                   LatGPM2[LatGPM2[,1]==NearGP(i,cGPM[,c(1,2)])[[2]],2] ])

MyGridGPM <- matrix(MyGridGPM[-1],dim(MyGrid)[1])
colnames(MyGridGPM) <- tGPM
rownames(MyGridGPM) <- 1:dim(MyGrid)[1]
dim(MyGridGPM)

# D) ERA5 DATA GRIDPOINTS -----------------------
cERA5 <- coordinates(pERA5)
cERA5 <- cbind(cERA5,1:dim(cERA5)[1])

MyGridERA5 <- 0
for(tt in 1:dim(pERA5)[2])
  for(i in 1:dim(MyGrid)[1])
    MyGridERA5 <- c(MyGridERA5, pERA5[[tt]][ 
                    cERA5[(cERA5[,1]==NearGP(i,cERA5[,c(1,2)])[[1]])&
                          (cERA5[,2]==NearGP(i,cERA5[,c(1,2)])[[2]]),3] ])

MyGridERA5 <- matrix(MyGridERA5[-1],dim(MyGrid)[1])
colnames(MyGridERA5) <- tERA5
rownames(MyGridERA5) <- 1:dim(MyGrid)[1]
dim(MyGridERA5)
MyGridERA5[1:10,1:5]


# MYGRID PRECIPITATION PLOT ====================================================

GridPLOT <- function(DATASET,TIME) {
  plot(c(0,0), type = "n", asp = 1.2, cex.main = 0.75,
    xlim = c(min(MyGrid[,1]),max(MyGrid[,1])), ylim = c(min(MyGrid[,2]),max(MyGrid[,2])),
    xlab = "Longitude (°E)", ylab = "Latitude (°N)",
    main = paste0(substitute(DATASET)," Precipitation\n",
               "La Spezia (Italy)\n",asDATE(colnames(DATASET)[TIME]) ),
    panel.first = c(ItalyBG(rgb(.95,.95,.95),"gray"), grid())   )
  .filled.contour(unique(MyGrid[,1]),unique(MyGrid[,2]),
    t(matrix(DATASET[,TIME],length(unique(MyGrid[,2])))),
    levels = c(.1,2,4,8,16,1000)*ifelse(substitute(DATASET)=="MyGridSAT",.0001,1),
    col = c(rgb(0,0,1,.5),rgb(0,1,0,.5),rgb(1,1,0,.5),rgb(1,.65,0,.5),rgb(1,0,0,.5)))
  points(MyGrid, pch = 3)
  PrecLEG(ifelse(substitute(DATASET)=="MyGridSAT","x10-4 kg/m?s","mm/h")) }

GridPLOT(MyGridSAT, 84)
GridPLOT(MyGridPLV, 22)
GridPLOT(MyGridGPM, 43)
GridPLOT(MyGridERA5,22)


# Check MyGrid Precipitation Plot ----------------------------------------------
PlotSAT(84,9.6,10.2,43.9,44.35,ItalyBG)
PlotPLV(22,9.6,10.2,43.9,44.35,ItalyBG)
PlotGPM(43,9.6,10.2,43.9,44.35,ItalyBG)
PlotERA5(22,9.6,10.2,43.9,44.35,ItalyBG)



# EXPORT PROCESSED DATA ========================================================


write.table(cbind(MyGrid,MyGridSAT),
            "MyGridSAT_LaSpezia.txt",
            sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(cbind(MyGrid,MyGridPLV),
            "MyGridPLV_LaSpezia.txt",
            sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(cbind(MyGrid,MyGridGPM),
            "MyGridGPM_LaSpezia.txt",
            sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(cbind(MyGrid,MyGridERA5),
            "MyGridERA5_LaSpezia.txt",
            sep = "\t", row.names=FALSE, col.names=TRUE)


# Save plot as PNG for GPM time serie (also accepts PDF, JPG, etc...) ----------
outDIR <- paste0("GPM_",tGPM[1],"_",tGPM[length(tGPM)])
dir.create(outDIR)
for(i in 1:length(tGPM)) {
  png(file=paste0(outDIR,"/GPM_",tGPM[i],".png"))
  PlotGPM(i,6,20,35,47,GlobalBG)
  dev.off()
}