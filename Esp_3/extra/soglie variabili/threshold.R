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
  "ggsci"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

statsFun <- function(Var, Ref) {
  avgR <- mean(Ref, na.rm = TRUE)
  SMP <- length(na.omit(Var))
  COR <- cor(Var, Ref, use = "pairwise.complete.obs")
  ME <- mean(Var - Ref, na.rm = TRUE)
  MAE <- mean(abs(Var - Ref), na.rm = TRUE)
  RMSE <- sqrt(mean((Var - Ref)^2, na.rm = TRUE))
  CV <- RMSE / avgR
  nME <- ME / avgR
  nMAE <- MAE / avgR
  MAX <- max(Var, na.rm = TRUE)

  table <- data.frame(THR, SMP, avgR, nME, nMAE, CV, COR, ME, MAE, RMSE, MAX)
  return(round(table, 2))
}

statsFun_CAT <- function(Var, Ref) {
  cat <- data.frame(
    Var = Var == 1,
    Ref = Ref == 1
  )

  xtab <- xtabs(~ Ref + Var, data = cat, na.action = "na.omit")

  hits <- as.double(xtab["TRUE", "TRUE"]) # [2,2]#
  crej <- as.double(xtab["FALSE", "FALSE"]) # [1,1]# #correct rejection
  miss <- as.double(xtab["TRUE", "FALSE"]) # [2,1]#
  fala <- as.double(xtab["FALSE", "TRUE"]) # [1,2]#

  FAR <- fala / (hits + fala)
  BIAS <- (hits + fala) / (hits + miss)
  POD <- hits / (hits + miss)
  CSI <- hits / (hits + miss + fala)
  SMP <- sum(xtab)
  ACC <- (hits + crej) / (hits + fala + miss + crej)
  POFD <- fala / (fala + crej)
  rand <- (hits + miss) * (hits + fala) / SMP
  ETS <- (hits - rand) / (hits + miss + fala - rand)

  table <- data.frame(THR, SMP, hits, crej, miss, fala, rand, POD, FAR, BIAS, ETS, CSI, ACC, POFD)
  return(round(table, 2))
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

figfolder <- "figure/"
resfolder <- "risultati/"

# Create the project directories, if missing
if (!dir.exists(figfolder)) {
  dir.create(figfolder)
}
if (!dir.exists(resfolder)) {
  dir.create(resfolder)
}

# Define time interval (between 1 October and 30 November 2016)
sel_START <- "201610120000"
sel_END <- "201610162300"

# Select properties of the common grid ("MyGrid"):
MyGrid_res <- 0.1
MyGrid_extLON <- c(0, 25)
MyGrid_extLAT <- c(35, 50)

# defining mask
GlobalBG <- vect(paste0(datafolder, "World_Countries_(Generalized)"))
ItalyBG <- GlobalBG[112]

# Set the unic conversion coefficients
CorrectionFactor_HSAF <- 3600 # kg/m2*s --> mm/h
CorrectionFactor_ERA5 <- 1000 # m/h --> mm/h

THRvector <- c(0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0) # mm/h

# initializing empty dataframes
emptydf     <- data.frame(THR = double(), SMP = double(), avgR = double(), 
                          nME = double(), nMAE = double(), CV = double(), 
                          COR = double(), ME = double(), MAE = double(), 
                          RMSE = double(), MAX = double())

emptydf_CAT <- data.frame(THR = double(), SMP = double(), hits = double(), 
                          crej = double(), miss = double(), fala = double(), 
                          rand = double(),  POD = double(), FAR = double(), 
                          BIAS = double(), ETS = double(), CSI = double(), 
                          ACC = double(), POFD = double()) 

ERA5stats <- emptydf
GPMstats  <- emptydf
HSAFstats <- emptydf
ERA5statsCAT <- emptydf_CAT
GPMstatsCAT  <- emptydf_CAT
HSAFstatsCAT <- emptydf_CAT

ERA5 <- rast(paste0(resfolder, "ERA5_final.tif"))
RG   <- rast(paste0(resfolder, "RG_final.tif"))
GPM  <- rast(paste0(resfolder, "GPM_final.tif"))
HSAF <- rast(paste0(resfolder, "HSAF_final.tif"))

names(RG)[which(!names(RG) %in% names(HSAF))]
new.times <- names(RG)[which(names(RG) %in% names(HSAF))]

RG   <- RG[[new.times]]
ERA5 <- ERA5[[new.times]]
GPM  <- GPM[[new.times]]

for (THR in THRvector) {
  # CREAZIONE DI UNICO DATASET =====================================
  ERA5.mask <- mask(ERA5 >= THR, ItalyBG)
  GPM.mask <- mask(GPM >= THR, ItalyBG)
  RG.mask <- mask(RG >= THR, ItalyBG)
  HSAF.mask <- mask(!is.na(HSAF), ItalyBG)

  # con corrispondenza 1:1 fra elementi, ma perdendo informazione temporale e spaziale. Li rendo vettori (1-dim!)
  data <- tibble(
    RainGauges = as.vector(values(clamp(RG, THR, values = FALSE))),
    ERA5 = as.vector(values(clamp(ERA5, THR, values = FALSE))),
    GPM = as.vector(values(clamp(GPM, THR, values = FALSE))),
    HSAF = as.vector(values(clamp(HSAF, THR, values = FALSE)))
  ) %>% na.omit()

  data_CAT <- tibble(
    RainGauges = as.vector(values(RG.mask)),
    ERA5 = as.vector(values(ERA5.mask)),
    GPM = as.vector(values(GPM.mask)),
    HSAF = as.vector(values(HSAF.mask))
  ) %>% na.omit()

  # STATISTICA SUI DATI ======================================
  # indicatori continui
  if (length(data$RainGauges) != 0) {
    if (length(data$ERA5) != 0) {
      ERA5stats <- rbind(ERA5stats, statsFun(data$ERA5, data$RainGauges))
    }
    if (length(data$GPM) != 0) {
      GPMstats <- rbind(GPMstats, statsFun(data$GPM, data$RainGauges))
    }
    if (length(data$HSAF) != 0) {
      HSAFstats <- rbind(HSAFstats, statsFun(data$HSAF, data$RainGauges))
    }
  }
  print(paste0("done ", THR))

  # indicatori discreti
  if (length(data_CAT$RainGauges) != 0) {
    if (length(data_CAT$ERA5) != 0) {
      ERA5statsCAT <- rbind(ERA5statsCAT, statsFun_CAT(data_CAT$ERA5, data_CAT$RainGauges))
    }
    if (length(data_CAT$GPM) != 0) {
      GPMstatsCAT <- rbind(GPMstatsCAT, statsFun_CAT(data_CAT$GPM, data_CAT$RainGauges))
    }
    if (length(data_CAT$HSAF) != 0) {
      HSAFstatsCAT <- rbind(HSAFstatsCAT, statsFun_CAT(data_CAT$HSAF, data_CAT$RainGauges))
    }
  }
  print(paste0("doneCAT ", THR))
}

# exporting to text files
# write.table(ERA5stats, paste0(resfolder, "ERA5stats.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
# write.table(GPMstats, paste0(resfolder, "GPMstats.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
# write.table(HSAFstats, paste0(resfolder, "HSAFstats.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
# write.table(ERA5statsCAT, paste0(resfolder, "ERA5statsCAT.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
# write.table(GPMstatsCAT, paste0(resfolder, "GPMstatsCAT.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
# write.table(HSAFstatsCAT, paste0(resfolder, "HSAFstatsCAT.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

# rearranging data in two dataframes for easier plots production
ERA5stats$SYS <- rep("ERA5", length(ERA5stats$THR))
GPMstats$SYS <- rep("GPM", length(GPMstats$THR))
HSAFstats$SYS <- rep("HSAF", length(HSAFstats$THR))
stats <- rbind(ERA5stats, GPMstats, HSAFstats)

ERA5statsCAT$SYS <- rep("ERA5", length(ERA5statsCAT$THR))
GPMstatsCAT$SYS <- rep("GPM", length(GPMstatsCAT$THR))
HSAFstatsCAT$SYS <- rep("HSAF", length(HSAFstatsCAT$THR))
statsCAT <- rbind(ERA5statsCAT, GPMstatsCAT, HSAFstatsCAT)

save(stats, statsCAT, file = paste0(resfolder, "thresholdBIS.RData"))


# plots

ggplot(stats, aes(THR, COR, color = SYS)) + geom_line() + theme_bw()
ggplot(stats, aes(THR, ME, color = SYS)) + geom_line() + theme_bw()
ggplot(stats, aes(THR, MAE, color = SYS)) + geom_line() + theme_bw()

ggplot(statsCAT, aes(THR, BIAS, color = SYS)) + geom_line() + theme_bw()
ggplot(statsCAT, aes(THR, CSI, color = SYS)) + geom_line() + theme_bw()
ggplot(statsCAT, aes(THR, ETS, color = SYS)) + geom_line() + theme_bw()


