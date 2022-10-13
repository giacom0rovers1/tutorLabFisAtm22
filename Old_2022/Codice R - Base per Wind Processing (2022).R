# ===========================================================
# = LABORATORY OF DATA PROCESSING AND EVALUATION            =
# = Measured Wind Speed and Wind Direction      v2020.10.05 =
# ===========================================================

# Clean variables previously defined
rm(list=ls())

# READING DATA ==============================================
# Define the working directory
setwd("C:/Users/qucci/Desktop/RTD_2017/Corso Lab/Corso 2021_22")

# Define the input filenames and assign them to a variable
StationDataFile <- "Dati Meteo - Stazione A.txt"

# Read data as a table and assigning to the variable called 'StationData'
StationData  <- read.table(file(StationDataFile),
                           sep       = "",              # separate columns by space
                           na.string = NA,              # tag empty values as "NA"
                           as.is     = TRUE,            # convert character variables to factors
                           header    = FALSE,           # ignore the header
                           skip      = 1,               # skip the first row
                           col.names = c("Date","Time","Interval","Press","inTemp","outTemp","inHum","outHum","wndSpeed","wndDir","windGust","wGstDir","rainRate","rain","dewpoint","wdChill","heatIndx","ET","totRad","UV")
                           )

#help(command)                     # Help about "command"
#ls()                              # Print the current active variables
#rm(var)                           # Clear the variable var
#A:B                               # Create a sequence from A to B
#seq_len(num)                      # Create a sequence of length num
#paste0("String",seq_len(num))     # Create a sequence of length num with a string
#substr(String,A,B)                # Take a piece of a string from the position A to B

# Print the dimension of a matrix
dim(StationData)

# Print only the first 6 rows and 10 columns of a matrix
StationData[1:6,1:10]


# DATE FORMAT ===============================================
# Define the first and second columns as a String and split the date into 5 columns
Temp <- (t(rbind(as.numeric(substr(StationData[,1],1,4)),  # year
                 as.numeric(substr(StationData[,1],5,6)),  # month
                 as.numeric(substr(StationData[,1],7,8)),  # day
                 floor(StationData[,2]),                   # hour
                 round(100*(StationData[,2] %%1))          # minutes
                 )))

# Define names to rows and columns
colnames(Temp) <- c("Year","Month","Day","Hour","Minutes")
rownames(Temp) <- 1:dim(Temp)[1]

# Replace the new date format (2 cols > 5 cols) and keep only the WindRose variables
StationData <- cbind(Temp,StationData[,9:10])

# Record the length of data
LenSD <- dim(StationData)[1]
StationData[1:6,]
LenSD


# DATA SORTING ==============================================
# Function to convert time to absolute time (days after 01/01/1970, R standard)
AbsTIME <- function(i) round(as.numeric(as.Date(paste0(StationData[i,1],
                       ifelse(StationData[i,2]<10,"0",""),StationData[i,2],
                       ifelse(StationData[i,3]<10,"0",""),StationData[i,3]),
                       "%Y%m%d"))+(StationData[i,4]/24)+((StationData[i,5]/60)*(1/24)),digits=2)
# Total data period in days
AbsTIME(dim(StationData)[1])-AbsTIME(1)

# Calculate the absolute time and add to StationData
Temp           <- matrix(AbsTIME(1:dim(StationData)[1]),dim(StationData)[1],1)
colnames(Temp) <- c("AbsTime")
StationData    <- cbind(Temp,StationData[,1:7])

StationData[1:6,]

# Sort StationData by the absolute time
StationData<-StationData[order(StationData$"AbsTime"),]
rownames(StationData) <- 1:dim(StationData)[1]

StationData[1:6,]
dim(StationData)

# DATA FILTERING ============================================

# Needs DPLYR Library
#install.packages("dplyr")
library(dplyr)
# KEEP ONLY HOURLY RECORDS ----------------------------------
# Create a variable to storage the excluded values
ExcTIME <- dplyr::filter(StationData,(StationData$"Minutes" > 0))

StationData <- dplyr::filter(StationData,(StationData$"Minutes" == 0))

# Upload the length of data
LenSD <- c(LenSD,dim(StationData)[1])
StationData[1:6,]
ExcTIME[1:6,]
LenSD


# FILTER RECORDS WITH WS<0 ----------------------------------
# Create a variable to storage the excluded values
ExcWSneg <- dplyr::filter(StationData,(StationData$"wndSpeed" < 0))

StationData <- dplyr::filter(StationData,(StationData$"wndSpeed" >= 0))

# Upload the length of data
LenSD <- c(LenSD,dim(StationData)[1])
StationData[1:6,]
ExcWSneg [1:6,]
LenSD


# FILTER RECORDS WITH WS CONSTANT FOR 3h --------------------
ExcWS3h <- 0     # Create a variable to storage the excluded values

# Find records under the condition
for (i in 3:dim(StationData)[1])
   if (((StationData[i,1]-StationData[i-2,1])<0.1)&    # 1h ~ 0.04 AbsTime
       (abs(StationData[i,7]-StationData[i-2,7])<0.1)&
       (abs(StationData[i,7]-StationData[i-1,7])<0.1)&
       (StationData[i,7]>0)  # keep calm winds
       )
      ifelse(length(ExcWS3h)==1,
         ExcWS3h <- c(i-2,i-1,i),
            ExcWS3h <- c(ExcWS3h,i-2,i-1,i))

# Exclude Duplicates (i)
ExcWS3h <- unique(ExcWS3h)

# Exclude conditioned records from StationData
Temp        <- ExcWS3h
ExcWS3h     <- StationData[Temp,]
StationData <- StationData[-Temp,]

# Upload the length of data
rownames(StationData) <- 1:dim(StationData)[1]
rownames(ExcWS3h)     <- 1:dim(ExcWS3h)[1]
LenSD                 <- c(LenSD,dim(StationData)[1])
StationData[1:6,]
ExcWS3h[1:6,]
LenSD

# EXPORT PROCESSED DATA =====================================
# Set the output filenames based on the input filename
StationDataFile <- strsplit(StationDataFile,".txt")

write.table(StationData, paste0(StationDataFile,"_proc.txt")    , sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(ExcTIME    , paste0(StationDataFile,"_ExcTIME.txt") , sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(ExcWSneg   , paste0(StationDataFile,"_ExcWSneg.txt"), sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(ExcWS3h    , paste0(StationDataFile,"_ExcWS3h.txt") , sep = "\t", row.names=FALSE, col.names=FALSE)


# DATA CONSISTENCY ===========================================
# Total of Hourly Data Expected by Year (Function)
yexp.data <- function(year) 24*(365 + ifelse(((year %% 400)/4) == floor((year %% 400)/4),1,0))
# Total of Hourly Data Expected by Month (Function)
mexp.data <- function(month) 24*(ifelse(month==4|month==6|month==9|month==11,30,
                                 ifelse(month==1|month==3|month==5|month==7|month==8|month==10|month==12,31,
                                 28 )))

# Print the summary of StationData
summary(StationData)

# Print the Data Consistency Table
Temp     <- c(sum(yexp.data(unique(StationData$"Year"))),
              LenSD[1], min(LenSD), max(LenSD)-min(LenSD),
              LenSD[1]-LenSD[2], LenSD[2]-LenSD[3], LenSD[3]-LenSD[4])
DataCONS <- matrix(as.numeric(c(Temp,round(Temp*100/sum(yexp.data(unique(StationData$"Year"))),digits=2)
                   )),7,2)
colnames(DataCONS) <- c("Station A","(%)")
rownames(DataCONS) <- c("Expected","Available","Valid","Filtered"," > by Time","  > by WSneg","   > by WS3h")
DataCONS

# DISPLAY THE DATA CONSISTENCY BY MONTH AS A BARPLOT ---------
# Compute the percentiles of valid data for each month
monPERC <- function(month) 100*dim(dplyr::filter(StationData,(StationData$"Month" == month)))[1]/
                                  (mexp.data(month)*length(unique(StationData$"Year"))+
                                   sum(ifelse(month == 2 & ((unique(StationData$"Year") %% 400)/4) ==
                                   floor((unique(StationData$"Year") %% 400)/4),1,0))   )

for (i in 1:12) ifelse(i==1,iPLOT<-monPERC(i),iPLOT[i]<-monPERC(i))

# Create a function for graphics as barplot
PLOT <- function(data) barplot(data, 
                       xlab = "Months", ylab = "Valid Data (%)", ylim = c(0,100), col = "darkblue",
                       main = paste0("Station A (",
                       ifelse(length(unique(StationData$"Year")) == 1, unique(StationData$"Year"),
                       paste0(min(unique(StationData$"Year")),"-",max(unique(StationData$"Year"))) ),
                       ")"),
                       names.arg = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")   )
PLOT(iPLOT)

# Save plot as PNG (also accepts PDF, JPG, etc...)
png(file=paste0(StationDataFile,"_DataConsistency.png"))
PLOT(iPLOT)
dev.off()


# HISTOGRAM OF WIND SPEED =====================================
HistWS <- function(data) hist(data, main = "Histogram of Wind Speed (Station A)",
                         xlab = "Wind Speed (m/s)", ylab = "Frequency",
                         density = 10, angle = 45, col = "#0066CC", labels = TRUE)

HistWS(StationData$"wndSpeed")

# Save plot as PNG (also accepts PDF, JPG, etc...)
png(file=paste0(StationDataFile,"_HistWS.png"))
HistWS(StationData$"wndSpeed")
dev.off()


# WIND ROSE ==================================================
WindROSE <- c(1:(8*5))  # 8 Wind Directions x 5 Wind Categories
wsCAT    <- c(2,6,12,19,30,500)/3.6
ii       <- 1           # WindROSE vector control

# NORTH WINDS LOOPING ----------------------------------------
for(ws in 2:6) {
   WindROSE[ii] <- dim(dplyr::filter(StationData,(
                      ((StationData$"wndDir" >= 360-22.5)|(StationData$"wndDir" < 0+22.5))&
                      (StationData$"wndSpeed" >= wsCAT[ws-1])&(StationData$"wndSpeed" < wsCAT[ws]))))[1]
   ii <- ii + 1 }

# OTHER DIRECTION WINDS LOOPING ------------------------------
for(wd in seq(45,315,45)) { for(ws in 2:6) {
   WindROSE[ii] <- dim(dplyr::filter(StationData,(
                      ((StationData$"wndDir" >= wd-22.5)&(StationData$"wndDir" < wd+22.5))&
                      (StationData$"wndSpeed" >= wsCAT[ws-1])&(StationData$"wndSpeed" < wsCAT[ws]))))[1]
      ii <- ii + 1 }}

# MATRIX BUILD -----------------------------------------------
WindROSE           <- 100*matrix(WindROSE,5,8)/sum(yexp.data(unique(StationData$"Year")))#dim(StationData)[1]
colnames(WindROSE) <- c("N","NE","E","SE","S","SW","W","NW")
rownames(WindROSE) <- c("0.56-1.67","1.67-3.33","3.33-5.28","5.28-8.33",">8.33")
WindROSE

# Get the stacked barplot ------------------------------------
WindRosePLOT <- function(data) {barplot(data, 
                xlab = "Wind Directions", ylab = "Frequency (%)", ylim = c(0,12), 
                border = "white", space = 0.05, font.axis = 2,
                col = c("#99CCFF","#6699CC","#0066CC","#003366","#000033"),
                main = paste0("Station A (",
                ifelse(length(unique(StationData$"Year")) == 1, unique(StationData$"Year"),
                paste0(min(unique(StationData$"Year")),"-",max(unique(StationData$"Year"))) ),")"),
                sub = paste0("Calm ",round(100*dim(dplyr::filter(StationData,StationData$"wndSpeed" < 2/3.6))[1]/
                                     sum(yexp.data(unique(StationData$"Year"))),digits=2),
                             "%     Missing ",round(100*(1-dim(StationData)[1]/sum(yexp.data(unique(StationData$"Year")))),
                digits=2),"%"))
                legend("topright",c(paste0("0.56 - 1.67 (",round(sum(WindROSE[1,]),digits=2),"%)"),
                            paste0("1.67 - 3.33 (",round(sum(WindROSE[2,]),digits=2),"%)"),
                            paste0("3.33 - 5.28 (",round(sum(WindROSE[3,]),digits=2),"%)"),
                            paste0("5.28 - 8.33 (",round(sum(WindROSE[4,]),digits=2),"%)"),
                            paste0("> 8.33 ("     ,round(sum(WindROSE[5,]),digits=2),"%)")),
                        pch   = 15, border = NA,
                        col   = c("#99CCFF","#6699CC","#0066CC","#003366","#000033"),
                        title = "Wind Speed (m/s)")}
WindRosePLOT(WindROSE)

# Save plot as PNG (also accepts PDF, JPG, etc...)
png(file=paste0(StationDataFile,"_WindRose.png"))
WindRosePLOT(WindROSE)
dev.off()


# EXPORT TO LAKE FORMAT ======================================

toLAKES <- cbind(matrix(rep(99999,dim(StationData)[1]),dim(StationData)[1],1), # Station ID (flag)
                 StationData[,2:5],                                            # Date (year, month, day, hour)
                 round(StationData[,8],digits=0),                              # Rounded Wind Direction (?)
                 round(1.94384*StationData[,7],digits=0)                       # Rounded Wind Speed (knot = 1.94384 x m/s)
                 )
colnames(toLAKES) <- c("StID","YYYY","MM","DD","HH","Dir","WS")
rownames(toLAKES) <- c(1:dim(StationData)[1])

StationData[1:6,]
toLAKES[1:6,]



#this section added to append the "LAKES FORMAT" string at the beginning of the file
fileConn<-file(paste0(StationDataFile,"_LAKES.txt"))
writeLines(c("LAKES FORMAT"), fileConn)
close(fileConn)

write.table(toLAKES, paste0(StationDataFile,"_LAKES.txt"),
            sep = " ", row.names=FALSE, col.names=FALSE, append =TRUE)

# STATISTICAL INDEXES ========================================
# Create a matrix with two columns containing the data to be compared
StatIND <- cbind(StationData[,7], runif(dim(StationData)[1], 0.5, 2)*StationData[,7])

# Bias
mean(StatIND[,1]-StatIND[,2])

# Root Mean Square Error (RMSE)
sqrt(sum((StatIND[,1]-StatIND[,2])^2)/dim(StatIND)[1])

# Pearson Correlation (R)
cor(StatIND)[1,2]
mean( (StatIND[,1]-mean(StatIND[,1]))*(StatIND[,2]-mean(StatIND[,2])) )/
                     (sd(StatIND[,1])*sd(StatIND[,2]))



# SCATTER PLOT AND LINEAR REGRESSION =========================
# Get the linear regression formula as text
TextLM <- strsplit(as.character(lm(StationData$"wndSpeed"~StationData$"Hour")[1])," ")
TextLM <- paste0("y = ",round(as.numeric(gsub(",","",TextLM[[1]][3])),digits=2),
                 " + ",round(as.numeric(gsub(")","",TextLM[[1]][6])),digits=2),"x")

# Save plot as PNG (also accepts PDF, JPG, etc...)
png(file=paste0(StationDataFile,"_ScatterPlot.png"))
plot(StationData$"wndSpeed"~StationData$"Hour", main = "Scatter Plot (Station A)",
     xlab = "Hour", ylab = "Wind Speed (m/s)", pch = 20, col = "gray")
abline(lm(StationData$"wndSpeed"~StationData$"Hour"), col = "red")
text(4,6, TextLM, col = "red")
dev.off()


# HOURLY MEAN WIND ============================================

WindHM <- function(hh) c(
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >= 360-22.5)| (StationData$"wndDir" <   0+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >=  45-22.5)& (StationData$"wndDir" <  45+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >=  90-22.5)& (StationData$"wndDir" <  90+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >= 135-22.5)& (StationData$"wndDir" < 135+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >= 180-22.5)& (StationData$"wndDir" < 180+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >= 225-22.5)& (StationData$"wndDir" < 225+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >= 270-22.5)& (StationData$"wndDir" < 270+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,(((StationData$"wndDir"   >= 315-22.5)& (StationData$"wndDir" < 315+22.5))&
                                (StationData$"wndSpeed" >=    2/3.6)& (StationData$"Hour" == hh))))[1],
dim(dplyr::filter(StationData,  (StationData$"wndSpeed" <     2/3.6)& (StationData$"Hour" == hh)))[1])

iWindHM <- matrix(c(1:(9*24)),9,24)
rownames(iWindHM) <- c("N","NE","E","SE","S","SW","W","NW","Calm")
colnames(iWindHM) <- c(0:23)
for(h in 0:23) {iWindHM[,h+1] <- 2*WindHM(h)/sum(WindHM(h))}  # x2 to normalize the scale with wind speed
iWindHM           <- matrix(iWindHM,9,24)

WindHM2  <- function(hh) mean(dplyr::filter(StationData,((StationData$"wndSpeed" >= 0*2/3.6)&(StationData$"Hour" == hh)))[,7])
iWindHM2 <- matrix(c(1:24),24,1)
for(h in 0:23) {iWindHM2[h+1,1] <- WindHM2(h)}

WindHM3  <- function(hh) sd(dplyr::filter(StationData,((StationData$"wndSpeed" >= 0*2/3.6)&(StationData$"Hour" == hh)))[,7])
iWindHM3 <- matrix(c(1:24),24,1)
for(h in 0:23) {iWindHM3[h+1,1] <- WindHM3(h)}

Colors <- c(rgb(219/255,236/255,246/255),rgb(201/255,228/255,240/255),rgb(184/255,218/255,237/255),
            rgb(165/255,209/255,232/255),rgb(198/255,227/255,219/255),rgb(209/255,233/255,227/255),
            rgb(221/255,238/255,233/255),rgb(232/255,244/255,241/255),rgb(250/255,250/255,250/255) )

WindHMPLOT  <- function(iplot) {
WindHMPLOTo <-  barplot(iplot, border = NA, space = 0, xlab = "Time (hours)",
                        ylab = "Wind Speed (m/s) - Mean (solid) + StDev (dashed)",
                        main = "Hourly Wind Speed and Wind Direction (Station A)",
                        names=c("00","","","03","","","06","","","09","","",12,"","",15,"","",18,"","",21,"",""),
                        col = Colors, ylim = c(0,2.2) )
abline(v = seq(3,21,3), col = "gray")
lines (x = WindHMPLOTo, y = iWindHM2)
points(x = WindHMPLOTo, y = iWindHM2, pch = 19)
lines (x = WindHMPLOTo, y = iWindHM3, lty = 2)
legend(.75, 2.2 , c("N","NE","E","SE","S","SW","W","NW","Calm"), pch = 15, border = NA, horiz = TRUE,
                title = "Frequency (%) by Wind Direction", col = Colors, bg = "white", cex = 0.8 ,pt.cex = 1.75)
}

# Save plot as PNG (also accepts PDF, JPG, etc...)
png(file=paste0(StationDataFile,"_DailyWind.png"))
WindHMPLOT(iWindHM)
dev.off()
dev.off()