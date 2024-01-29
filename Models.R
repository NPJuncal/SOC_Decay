setwd("C:/Users/npjun/Dropbox/Seagrasses/Degradacion anaerobia_OSCAR/SOC_Decay/Data")

library(rbacon)
library(ggplot2)
library(dplyr)


### chronological models estimation

File<-"Dates.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ";",
                  dec = ".",
                  fileEncoding="latin1")
Dates<-as.data.frame(Dates)


X<-split(Dates, Dates$Core)

models <- data.frame(Core=character(),
                 Depth=numeric(), 
                 Acc.Mass=numeric(),
                 Age=character(),
                 Age.Pb=character(),
                 Age.C=character()) 



for(i in 1:length(X)) {
  
  
  Data<-as.data.frame(X[i])
  colnames(Data)<-colnames(Dates)
  
  
  Bacon_D <- as.data.frame(matrix(NA, nrow = nrow(Data), ncol = 7))
  colnames(Bacon_D)<-c("labID", "age", "error","depth","cc","dR","dSTD")
  SamY<-as.data.frame(matrix(NA, ncol = 7))
  colnames(SamY)<-c("labID", "age", "error","depth","cc","dR","dSTD")
  
  Bacon_D$labID<-paste(Data$Core.ID,Data$Depth)
  Bacon_D$age<-Data$Raw.BP
  Bacon_D$error<-Data$Error
  Bacon_D$depth<-Data$Acc.Mass
  Bacon_D$cc<-Data$Bacon.cc
  Bacon_D$dR<-Data$Bacon.dr
  Bacon_D$dSTD<-Data$Bacon.dSTD
  
  SamY$labID<-paste(Data$Core.ID[1],"0")
  SamY$age<-1950-(Data$Sampling.Y[1])
  SamY$error<-"1"
  SamY$depth<-"0"
  SamY$cc<-"0"
  SamY$dR<-"0"
  SamY$dSTD<-"0"
  
  
  Bacon_D<-na.omit(Bacon_D)
  Bacon_Pb<-subset(Bacon_D, cc == 0) #to select rows with Pb values
  Bacon_14C<-subset(Bacon_D, cc == 2 | cc == 1) #to select rows with NO Pb values (in our data 14C values)
  Bacon_14C<-rbind(SamY, Bacon_14C) # to add the minimum age
  
  
  ###create folders and save bacon files ####
  
  dir.create("Bacon_runs")
  path <- "Bacon_runs"
  
  ###create data frame with depth, DBD, POC and bacon models ####
  
  temp <-data.frame(matrix(NA, nrow = nrow(Data), ncol = 6))
  colnames(temp)<-c("Core","Depth","Acc.Mass","Age", "Age.Pb", "Age.C")
  
  
  temp$Core<-Data$Core#as depth we chose the center of the sample
  temp$Depth<-Data$Depth
  temp$Acc.Mass<-Data$Acc.Mass
  
  ## if column Age of the data has values, this is a cronological model from 210Pb, we put it in Pb column and skip model estimation because we dont have the raw data
  if (all(is.na(Data$Age))==FALSE) {
    temp$Age.Pb<-Data$Age
    
  } else {
  
  
  #### estimation Pb-C models
  
  if("210Pb" %in% unique(Data$Method)  == TRUE & "14C" %in% unique(Data$Method)  == TRUE) 
 
  {newfolder<-Data$Core.ID[1]
  newpath<- file.path(path, newfolder)
  dir.create(newpath)
  write.csv(Bacon_D,file.path(newpath,paste0(newfolder,".csv")),row.names=FALSE, quote=F) # use paste0() no paste() to create file names o bacon will have trouble with the space between the file name and de .csv 
  
  #bacon with both
  bfile<-Data$Core.ID[1]
  Bacon(bfile,
        d.min=0,d.max=max(Data$Acc.Mass))
  
  #load chronological model file
  results<-list.files(path = file.path(path, Data$Core.ID[1]), full.names = TRUE)
  Ages.file<-results[grep("*_ages*", results)]
  crono1<-read.table(Ages.file, header = T, sep = "", dec = ".")
  
  # put data into temp data.frame
  for (j in 1:(length(temp$Acc.Mass))) {
    
    temp[j,"Age"]<-crono1[which.min(abs( crono1[,"depth"]- temp[j,"Acc.Mass"])),"mean"]
  }}
  
  
  #### estimation Pb models
  if("210Pb" %in% unique(Data$Method)  == TRUE) 
  
  {newfolder<-paste(Data$Core.ID[1],"Pb")
  newpath<- file.path(path, newfolder)
  dir.create(newpath)
  write.csv(Bacon_Pb,file.path(newpath,paste0(newfolder,".csv")),row.names=FALSE, quote=F)
  
  #bacon with Pb
  bfile<-paste(Data$Core.ID[1],"Pb")
  Bacon(bfile,
        d.min=0,d.max=max(Data$Acc.Mass))
  
  #load chronological model file
  results<-list.files(path = file.path(path, paste(Data$Core.ID[1],"Pb")), full.names = TRUE)
  Ages.file<-results[grep("*_ages*", results)]
  crono2<-read.table(Ages.file, header = T, sep = "", dec = ".") 
  
  # put data into temp data.frame
  for (j in 1:(length(temp$Acc.Mass))) {
    
    temp[j,"Age.Pb"]<-crono2[which.min(abs(crono2[,"depth"]- temp[j,"Acc.Mass"])),"mean"]
  }}
  
  
  #bacon with C
  if("14C" %in% unique(Data$Method)  == TRUE)
  {newfolder<-paste(Data$Core.ID[1],"C")
  newpath<- file.path(path, newfolder)
  dir.create(newpath)
  write.csv(Bacon_14C,file.path(newpath,paste0(newfolder,".csv")),row.names=FALSE, quote=F)
  
  #bacon with 14C
  bfile<-paste(Data$Core.ID[1],"C")
  Bacon(bfile,
        d.min=0,d.max=max(Data$Acc.Mass))
  
  #load chronological model file
  results<-list.files(path = file.path(path, paste(Data$Core.ID[1],"C")), full.names = TRUE)
  Ages.file<-results[grep("*_ages*", results)]
  crono3<-read.table(Ages.file, header = T, sep = "", dec = ".")
  
  #we fill the chronological model columns with the mean age at the Acc mass 
  
  for (j in 1:(length(temp$Acc.Mass))) {
    
    temp[j,"Age.C"]<-crono3[which.min(abs(crono3[,"depth"]- temp[j,"Acc.Mass"])),"mean"]
  }}}
  
  
  models <- rbind(models, temp) }

## from years BP to years from sampling
list<-colnames(models)

models$SY<-NA

for (i in 1:nrow(models)) {
  
  models[i,"SY"]<-Dates[which(Dates$Core==models[i,"Core"])[1],"Sampling.Y"]
}

models <- models %>% mutate (Age1=abs((1950-SY)-Age)) 
models <- models %>% mutate (Age.Pb1=abs((1950-SY)-Age.Pb)) 
models <- models %>% mutate (Age.C1=abs((1950-SY)-Age.C)) 

models<-models[,c(1:3,8:10)]
colnames(models)<-list


write.csv(models,"Acc_Mass-Age_Gr.csv",sep=";", dec=",")  


#### add ages from new salt marsh cores

File <- "Data/Acc_Mass-Age_F.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ",",
                  dec = ".")
Dates <- as.data.frame(Dates)

