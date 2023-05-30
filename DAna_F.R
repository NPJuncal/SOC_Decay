setwd("C:/Users/npjun/Dropbox/Seagrasses/Degradacion anaerobia_OSCAR/SOC_Decay")

library(ggplot2)
library(maps)
library(tidyr)
library(magrittr)
library(mosaic)
library(gridExtra)
library(nlstools)
library(dplyr)

############## Loading the dataset #####################

File <- "DAta/Cores.csv"

Cores <- read.csv(File,
                  header = T,
                  sep = ";",
                  dec = ".")
Cores <- as.data.frame(Cores)

#clean the rows without OC data
Cores$Corg_.gcm3 <- as.numeric(Cores$Corg_.gcm3)
Cores <- Cores[!is.na(Cores$Corg_.gcm3),]

### we create a folder to save the results

Folder = "Decay2023"
dir.create(Folder)


### We select only vegetated or un-vegetated (Bare) cores

unique(Cores[, 5])

B = filter(Cores, Cores$V.vs.B != "Bare")
#Bare <- subset(A, A[,5]=='Bare')
length(unique(B$Core))


B$Corg <- as.numeric(B$Corg)
B$Mud <- as.numeric(B$Mud)
B$DBD <- as.numeric(B$DBD)

###### Sampling sites Map #########


B$Long <- as.numeric(B$Long)
B$Lat <- as.numeric(B$Lat)


#load a world map
WM <- map_data("world")

B %>%
  ggplot() + ggtitle("Sampling sites") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
  coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
  scale_fill_manual(values = c("blue",  "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(path = Folder,
  filename =  "Sampling sites.jpg",
  units = "cm",
  width = 20,
  height = 10
)




##### Average and median OC content per core (full length and top 25 cm)##################

ADT <- B[, c("Core", "Ecosystem", "Max.Depth", "Corg", "Mud")]

X <- split(ADT, ADT$Core)

CM <- data.frame(
  ID = character(),
  Ecosystem = character(),
  Av_C = numeric(),
  SE_C = numeric(),
  M_C = numeric(),
  Av_Mud = numeric(),
  SE_Mud = numeric(),
  M_Mud = numeric(),
  Av_C_25 = numeric(),
  SE_C_25 = numeric(),
  Av_Mud_25 = numeric(),
  SE_Mud_25 = numeric(),
  M_Mud_25 = numeric()
)

std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

for (i in 1:length(X)) {
  CM[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(ADT)
  CM[i, 2] <- Data[1, which(colnames(Data) == "Ecosystem")]
  #Data<-na.omit(Data)
  CM[i, 3] <- mean(Data[, which(colnames(Data) == "Corg")], na.rm=TRUE)
  CM[i, 4] <- std(Data[, which(colnames(Data) == "Corg")])
  CM[i, 5] <- median(Data[, which(colnames(Data) == "Corg")], na.rm=TRUE)
  CM[i, 6] <- mean(Data[, which(colnames(Data) == "Mud")], na.rm=TRUE)
  CM[i, 7] <- std(Data[, which(colnames(Data) == "Mud")])
  CM[i, 8] <- median(Data[, which(colnames(Data) == "Mud")], na.rm=TRUE)
  Data25 <-
    Data %>% filter(Data[, which(colnames(Data) == "Max.Depth")] <= 25)
  CM[i, 9] <- mean(Data25[, which(colnames(Data25) == "Corg")], na.rm=TRUE)
  CM[i, 10] <- std(Data25[, which(colnames(Data25) == "Corg")])
  CM[i, 11] <- mean(Data25[, which(colnames(Data25) == "Mud")], na.rm=TRUE)
  CM[i, 12] <- std(Data25[, which(colnames(Data25) == "Mud")])
}


ggplot(CM, aes(Ecosystem, Av_Mud_25)) +
  geom_boxplot() +
  geom_jitter()

ggplot(CM, aes(Av_C, Av_C_25))+
  geom_point()
ggplot(CM, aes(Av_Mud, Av_Mud_25))+
  geom_point()

ggplot(CM, aes(Av_Mud, Av_C))+
  geom_point()


CM %>%  group_by(Ecosystem) %>%
  summarise_at(vars(Av_C_25), list(name = median))

pairwise.wilcox.test(CM$Av_C_25, CM$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)
pairwise.wilcox.test(CM$Av_Mud_25, CM$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

write.csv(CM,
          file.path(Folder, "AvMdC.csv"),
          sep = ";",
          dec = ".")


# Sediment accretion rate -------------------------------------------------

File <- "Data/Acc_Mass-Age.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ",",
                  dec = ".")
Dates <- as.data.frame(Dates)

#SAR 

SAR <- data.frame(
  ID = character(),
  SAR_Age = numeric(),
  SAR_Pb = numeric(),
  SAR_C = numeric())

X<-split(Dates, Dates$Core)

for (i in 1:length(X)) {
  SAR[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(Dates)
  
  #correlation depth-age to predict depth at 150 yr old
  DataAge <- Data[!is.na(Data$Age),]
  DataAge<- DataAge[c(1:(length(which(DataAge$Age <=150)))),]
  
  if (nrow(DataAge)>2) {SAR[i, 2] <- max(DataAge$Depth)/max(DataAge$Age)}
  
  DataPb <- Data[!is.na(Data$Age.Pb),]
  DataPb<- DataPb[c(1:(length(which(DataPb$Age.Pb <=150)))),]
  
  if (nrow(DataPb)>2) {SAR[i, 3] <- max(DataPb$Depth)/max(DataPb$Age.Pb)}
  
  DataC <- Data[!is.na(Data$Age.C),]
  DataC<- DataC[c(1:(length(which(DataC$Age.C <=150)))),]
  
  if (nrow(DataC)>2) {SAR[i, 4] <- max(DataC$Depth)/max(DataC$Age.C)}
  
}

SAR$SAR<-"NA"
SAR$SAR<-as.numeric(SAR$SAR)

for (i in 1:nrow(SAR)) {
  
  if (is.na(SAR[i,"SAR_Age"]) == FALSE) {SAR[i,"SAR"]<-SAR[i,"SAR_Age"]} 
  
  else {
    
    if (is.na(SAR[i,"SAR_Pb"]) == FALSE) {SAR[i,"SAR"]<-SAR[i,"SAR_Pb"]}
    
   else {if (is.na(SAR[i,"SAR_C"]) == FALSE) {SAR[i,"SAR"]<-SAR[i,"SAR_C"]}}}}



# Corg trends with depth --------------------------------------------------

#### spearman correlation between depth and Corg content, as data no normal nor expected to follow lineal correlation
#### return DT data.frame (Core ID, rho, p value, trend)

ADT <- B[, c("Core", "Ecosystem", "Max.Depth", "Corg")]

X <- split(ADT, ADT$Core)

DT <- data.frame(
  ID = character(),
  Ecosystem = character(),
  C_rho = numeric(),
  C_p = numeric(),
  C_Gr = character()
)



for (i in 1:length(X)) {
  DT[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(ADT)
  DT[i, 2] <- Data[1, which(colnames(Data) == "Ecosystem")]
  #Data<-na.omit(Data)
  cor <-
    cor.test(x = Data[, which(colnames(Data) == "Max.Depth")], y = Data[, which(colnames(Data) ==
                                                                              "Corg")], method = "spearman")
  DT[i, 3] <- cor$estimate
  DT[i, 4] <- cor$p.value
  
}


##### from spearman correlations we discriminate 3 groups: NT, no trend with depth; DEC, decrease with depth; INC, increase with depth
#lower keys: IDs, capital keys: data frame with core data
# return data frames per group (NT,DEC and INC) and add trend tipe to DT

Nt <- DT$ID[DT$C_p > 0.5]
Trend <- DT[!(DT$C_p > 0.5), ]
Dec <- Trend$ID[Trend$C_rho < 0]
Inc <- Trend$ID[Trend$C_rho > 0]

NT <- B[is.element(B$Core, Nt), ]
DEC <- B[is.element(B$Core, Dec), ]
INC <- B[is.element(B$Core, Inc), ]

for (i in 1:nrow(DT)) {
  if (is.element(DT[i, 1], Nt)) {
    DT[i, 5] <- "NT"
  }
  else if (is.element(DT[i, 1], Dec)) {
    DT[i, 5] <- "DEC"
  }
  else {
    DT[i, 5] <- "INC"
  }
}


write.csv(DT,
          file.path(Folder, "DepthTrends.csv"),
          sep = ";",
          dec = ".")

##### Av Corg and trends

ggplot(CM, aes(DT$C_Gr, CM$Av_Mud)) +
  geom_boxplot() +
  geom_jitter()


#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

shapiro.test(CM$Av_Mud) #(>0.05 normal, <0.05 no normal)

## Student's t-test  if normally distributed, wilcox if not

pairwise.wilcox.test(CM$Av_Mud, DT$C_Gr,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

ggplot(CM, aes(DT$Ecosystem, CM$Av_Mud)) +
  geom_boxplot()

#### Count cores per group and get the percentages

NGr <- DT %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

DT %>% group_by(Ecosystem, C_Gr) %>% count()
NGrSg <- subset(DT, Ecosystem == "Seagrass") %>% group_by(C_Gr) %>% count()
NGrSg %>% mutate(proc = ((n * 100) / sum(NGrSg[, 2])))

NGrSm <-
  subset(DT, Ecosystem == "Salt Marsh") %>% group_by(C_Gr) %>% count()
NGrSm %>% mutate(proc = ((n * 100) / sum(NGrSm[, 2])))

NGrMg <- subset(DT, Ecosystem == "Mangrove") %>% group_by(C_Gr) %>% count()
NGrMg %>% mutate(proc = ((n * 100) / sum(NGrMg[, 2])))

### plot per grupos. Change size of jpg file when saving!!!!!!

ggplot(NT, aes(Max.Depth, Corg)) + xlab("Depth (cm)") + ylab("Organic carbon (g cm-3)") +
  geom_point(aes(Max.Depth, Corg)) +
  geom_line(aes(Max.Depth, Corg)) +
  facet_wrap( ~ Core, ncol = 5 , scales = "free") +
  coord_flip() +
  scale_x_reverse() +
  theme_light()

ggsave(
  path = Folder,
  filename = 'NTpor.jpg',
  width = 20,
  height = 200,
  units = 'cm',
  limitsize = FALSE
)

ggplot(DEC, aes(Max.Depth, Corg)) + xlab("Depth (cm)") + ylab("Organic carbon (g cm-3)") +
  geom_point(aes(Max.Depth, Corg)) +
  geom_line(aes(Max.Depth, Corg)) +
  facet_wrap( ~ Core, ncol = 5, scales = "free") +
  coord_flip() +
  scale_x_reverse() +
  theme_light()

ggsave(
  path = Folder,
  filename = 'DECpor.jpg',
  width = 20,
  height = 400,
  units = 'cm',
  limitsize = FALSE
)


ggplot(INC, aes(Max.Depth, Corg)) + xlab("Depth (cm)") + ylab("Organic carbon (g cm-3)") +
  geom_point(aes(Max.Depth, Corg)) +
  geom_line(aes(Max.Depth, Corg)) +
  facet_wrap( ~ Core, ncol = 5 , scales = "free") +
  coord_flip() +
  scale_x_reverse() +
  theme_light()

ggsave(
  path = Folder,
  filename = 'INCpor.jpg',
  width = 20,
  height = 200,
  units = 'cm',
  limitsize = FALSE
)


#### MDA ### 

CMS<-CM[,c(1,9,11)]
SARS<-SAR[,c(1,5)]
BS<- B[ !duplicated(B$Core), ]
names(BS)[names(BS) == 'Core'] <- 'ID'


SUM<-left_join(DT, CMS, by="ID")
SUM<-left_join(SUM, SARS, by="ID")
SUM<-left_join(BS, SUM, by="ID")

SUM<-SUM[,-c(15:19)]








###################################################################################
################ Correlation with time ############################################
###################################################################################

### chronological models estimation

File <- "Data/Acc_Mass-Age.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ",",
                  dec = ".")
Dates <- as.data.frame(Dates)

# Create a dataframe with those core with age-acc mass models


X <- split(B, B$Core)
C <- data.frame(matrix(ncol = 26, nrow = 0))
colnames(C) <- c(colnames(B), colnames(Dates[, 2:6]))

for (i in 1:length(X)) {
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(B)
  
  if (Data[1, 1] %in% Dates$Core == TRUE) {
    temp <- subset(Dates, Dates$Core == Data[1, 1])
    temp <- temp[c(1:nrow(Data)), ]
    Data <- cbind(Data, temp[, c(4:7)])
    
    C <- rbind(C, Data)
    
  } else
    next
}

length(unique(C$Core))


# df with core with any model

TAll = filter(C, !is.na(Age) | !is.na(Age.Pb) | !is.na(Age.C))
length(unique(TAll$Core))


#### Homogenize Age


TAll$FAge<-"NA"
TAll$FAge<-as.numeric(TAll$FAge)

for (i in 1:nrow(TAll)) {
  
  if (is.na(TAll[i,"Age"]) == FALSE) {TAll[i,"FAge"]<-TAll[i,"Age"]} 
  
  else {
    
    if (is.na(TAll[i,"Age.Pb"]) == FALSE) {TAll[i,"FAge"]<-TAll[i,"Age.Pb"]}
    
    else {if (is.na(TAll[i,"Age.C"]) == FALSE) {TAll[i,"FAge"]<-TAll[i,"Age.C"]}}}}


## mean depth at 150 years


X <- split(TAll, TAll$Core)

MaxDepth <- data.frame(ID = character(),
                       MaxDepth = numeric())

for (i in 1:length(X)) {
  MaxDepth[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  MaxDepth[i, 2] <- max(Data[, 15])
  
}

hist(MaxDepth$MaxDepth)
summary(MaxDepth$MaxDepth)


#### spearman correlation between time and Corg content
# for the last 150 years!

# function to generate different df depending on their tendency with time (DEC, INC and NT), generate figures per group 
# and return a dfs: TDT, TDEC, TINC and TNT. Needs a df with columns: c("Core", "Ecosystem","Min.Depth","Max.Depth", "FAge", "Corg"
# a max age (MA) and a name for the plots 

tendency<- function (df, MA, fn) {
  
  ADT <- df[, c("Core", "Ecosystem","Min.Depth","Max.Depth", "FAge", "Corg")]
  ADT<-subset(ADT, ADT$FAge < MA)
  
  X <- split(ADT, ADT$Core)
  
  TDT <- data.frame(
    ID = character(),
    Ecosystem = character(),
    C_rho = numeric(),
    C_p = numeric(),
    C_Gr = character()
  )  
  
  for (i in 1:length(X)) {
    TDT[i, 1] <- names(X[i])
    Data <- as.data.frame(X[i])
    colnames(Data) <- colnames(ADT)
    TDT[i, 2] <- Data[1, which(colnames(Data) == "Ecosystem")]
    #Data<-na.omit(Data)
    if (nrow(Data)<3) next
    cor <- cor.test(x = Data$FAge, y = Data$Corg, method = "spearman")
    TDT[i, 3] <- cor$estimate
    TDT[i, 4] <- cor$p.value
    
  }
  
  
  ##### from spearman correlations we discriminate 3 groups: NT, no trend with depth; DEC, decrease with depth; INC, increase with depth
  #lower keys: IDs, capital keys: data frame with core data
  # return data frames per group (NT,DEC and INC) and add trend tipe to DT
  
  
  Nt <- TDT$ID[TDT$C_p > 0.5]
  Trend <- TDT[!(TDT$C_p > 0.5), ]
  Dec <- Trend$ID[Trend$C_rho < 0]
  Inc <- Trend$ID[Trend$C_rho > 0]
  
  TNT <- df[is.element(df$Core, Nt), ]
  TDEC <- df[is.element(df$Core, Dec), ]
  TINC <- df[is.element(df$Core, Inc), ]
  
  for (i in 1:nrow(TDT)) {
    if (is.element(TDT[i, 1], Nt)) {
      TDT[i, 5] <- "NT"
    }
    else if (is.element(TDT[i, 1], Dec)) {
      TDT[i, 5] <- "DEC"
    }
    else {
      TDT[i, 5] <- "INC"
    }
  }
  
  #### Count cores per group and get the percentages
  
  library(dplyr)
  NGr <- TDT %>% group_by(C_Gr) %>% count()
  NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))
  
  # plot groups for visual check
  ggplot(TNT, aes(FAge, Corg)) + xlab("Time (years)") + ylab("Organic carbon (g cm-3)") +
    geom_point(aes(FAge, Corg)) +
    geom_line(aes(FAge, Corg)) +
    facet_wrap( ~ Core, ncol = 5 , scales = "free") +
    coord_flip() +
    scale_x_reverse() +
    theme_light()
  
  ggsave(
    path = Folder,
    filename = paste(fn,'_TNT.jpg'),
    width = 20,
    height = 30,
    units = 'cm'
  )
  
  ggplot(TDEC, aes(FAge, Corg)) + xlab("Time (years)") + ylab("Organic carbon (g cm-3)") +
    geom_point(aes(FAge, Corg)) +
    geom_line(aes(FAge, Corg)) +
    facet_wrap( ~ Core, ncol = 5, scales = "free") +
    coord_flip() +
    scale_x_reverse() +
    theme_light()
  
  ggsave(
    path = Folder,
    filename = paste(fn,'_TDEC.jpg'),
    width = 20,
    height = 50,
    units = 'cm',
    limitsize = FALSE
  )
  
  
  ggplot(TINC, aes(FAge, Corg)) + xlab("Time (years)") + ylab("Organic carbon (g cm-3)") +
    geom_point(aes(FAge, Corg)) +
    geom_line(aes(FAge, Corg)) +
    facet_wrap( ~ Core, ncol = 5 , scales = "free") +
    coord_flip() +
    scale_x_reverse() +
    theme_light()
  
  ggsave(
    path = Folder,
    filename = paste(fn,'_TINC.jpg'),
    width = 20,
    height = 30,
    units = 'cm'
  )
  
  return(list(TDT, TDEC, TINC, TNT))
}
  
prueba2<-tendency(TAll, MA=150, fn="prueba")









### groups and SAR ###

Prueba<-left_join(TDT,SAR, "ID")

ggplot(Prueba, aes(C_Gr, SAR))+
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))

ggplot(Prueba, aes(C_Gr, SAR))+
  geom_boxplot(aes(color=Ecosystem))+
  geom_jitter(aes(color=Ecosystem))


pairwise.wilcox.test(Prueba$SAR, Prueba$C_Gr,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


#################################
### NLM time-Accumulated Mass ###
#################################

#For those cores that decrease with time

DataA <-
  as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])


#Estimate organic carbon accumulated mass

estimate_h <- function(df = NULL) {
  
  # create individual data frames per each core
  
  df$Core <- factor(df$Core, levels=unique(df$Core))
  X<-split(df, df$Core)
  
  
  columns<-c("EMin","EMax","h")
  Fdf2 = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(Fdf2) = columns
  
  
  for(i in 1:length(X)) {
    
    Data<-as.data.frame(X[i])
    colnames(Data)<-colnames(df)
    
    #check if there is spaces between samples (e.g, first sample ends at 5 cm and next starts at 7)
    space<- c()
    
    for (j in 1:(nrow(Data)-1)) {
      
      # if there are no spaces between samples min and maximun depth of samples remain the same
      if (Data[j,which(colnames(Data)=="Max.Depth")] == Data[j+1,which(colnames(Data)=="Min.Depth")]) {
        space[j]<-FALSE} else {space[j]<-TRUE}}
    
    if (any(space==TRUE)) {
      # if there are spaces between samples it estimate the medium point between the maximum depth of the sample and the minimum
      #depth of the next sample and divide that distance between both samples
      Data <- cbind(Data, EMin=NA, EMax=NA)
      Data[1,"EMin"]<-0
      Data[nrow(Data),"EMax"]<-Data[nrow(Data),"Max.Depth"]
      for (j in 1:(nrow(Data)-1)) {
        if(space[j]==TRUE) {
          Data[j,"EMax"]<-Data[j,"Max.Depth"]+((Data[j+1,"Min.Depth"]-Data[j,"Max.Depth"])/2)
          Data[j+1,"EMin"]<-Data[j,"Max.Depth"]+((Data[j+1,"Min.Depth"]-Data[j,"Max.Depth"])/2)} else {
            Data[j,"EMax"]<-Data[j,"Max.Depth"]
            Data[j+1,"EMin"]<-Data[j+1,"Min.Depth"]}}
      
    }  else{
      Data <- cbind(Data, EMin=NA, EMax=NA)
      Data$EMin<-Data$Min.Depth
      Data$EMax<-Data$Max.Depth
      
    }
    
    Data <- cbind(Data, h=NA)
    
    #estimation of the thickness of the sample (h) from the new minimun and max depth of the sample
    
    Data<- Data |> dplyr::mutate (h = EMax-EMin)
    
    temp<-cbind(Data$EMin, Data$EMax, Data$h)
    colnames(temp)<-colnames(Fdf2)
    Fdf2<-rbind(Fdf2, temp)
    
  }
  Fdf<-cbind(df, Fdf2)
  
  return(Fdf)
}

DataA<-estimate_h(DataA)


#estimate carbon density and acc mass per sample

#organic carbon mass per sample

DataA<- DataA %>% mutate (OCg = DBD*(Corg/100)*h)

#Acc organic matter

DataAM<- DataA[0,]
DataAM[1,]=NA  # ad a temporary new row of NA values
DataAM[,'Corg.M'] = NA # adding new column, called for example 'new_column'
DataAM = DataAM[0,]

X<- split(DataA, DataA$Core)

for (i in 1:length(X)) {
  Data <- as.data.frame(X[i])
  colnames(Data)<-colnames(DataA)
  
  Data <- cbind(Data, Corg.M=NA)
  
  Data[1,"Corg.M"]<-Data[1,"OCg"]
  
  for (j in 2:nrow(Data)){
    Data[j,"Corg.M"]<-Data[j,"OCg"]+Data[j-1,"Corg.M"]
  }

  DataAM<-rbind(DataAM,Data)}

  
 #estimate model time acc mat 

DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])

OCModel<-function (df, nwpath) {
  
  SSS <- split(df, df$Core)
  
  
  Tfits <- data.frame(
    ID = character(),
    Ecosystem = character(),
    P = numeric(),
    k = numeric(),
    Max.Age = numeric()
  )
  
  cidr <- getwd()
  dir.create(file.path(cidr, nwpath), recursive = TRUE)  
  
  for (i in 1:length(SSS)) {
    Tfits[i, 1] <- names(SSS[i])
    Pr <- as.data.frame(SSS[i])
    colnames(Pr) <- list("Core", "Ecosystem", "FAge", "Corg", "Corg.M")
    Tfits[i, 2] <- Pr[1, which(colnames(Data) == "Ecosystem")]
    Tfits[i, 5] <- max(Pr$FAge)
    
    skip_to_next <- FALSE
    
    tryCatch(
      Exp1 <-
        nls(
          Corg.M ~ (p / k) * (1 - (exp(-k * FAge))),
          data = Pr,
          start = list(p = 0.01, k = 0.03)
        )
      ,
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    
    if (skip_to_next) {
      next
    }
    
    Func <-
      fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * FAge))),
               data = Pr,
               start = list(p = 0.01, k = 0.03))
    
    Coef <- coef(Exp1)
    
    finales <- as.list(coef(Exp1))
    Tfits[i, 3] <- as.numeric(finales[1])
    Tfits[i, 4] <- as.numeric(finales[2])
    
    
    fitY1 <- as.data.frame(c(1:max(Pr$FAge)))
    fitY1['new_col'] <- NA
    fitY1[, 2] <- Func(c(1:max(Pr$FAge)), Coef[1], Coef[2])
    colnames(fitY1) <- list("FAge", "Predict")
    
    p1 <-
      ggplot(Pr, aes(FAge, Corg)) + xlab("Age (years)") + ylab("Corg (g cm-3)") +
      geom_point() +
      geom_line() +
      coord_flip() +
      scale_x_reverse()
    
    p2 <-
      ggplot(fitY1, aes(FAge, Predict)) + xlab("Age (years)") + ylab("Corg Aumulated mass (g cm-2)") +
      geom_line(color = "blue") +
      geom_point(data = Pr, aes(FAge, Corg.M)) +
      coord_flip() +
      scale_x_reverse() +
      annotate(
        "text",
        x = max(Pr$FAge)/4,
        y = 0.01,
        label = "M=(p/k)*(1-exp(-k*Age)",
        hjust = "left"
      ) +
      annotate(
        "text",
        x = (max(Pr$FAge)/4)*2,
        y = 0.01,
        label = paste("p=", Coef[1]),
        hjust = "left"
      ) +
      annotate(
        "text",
        x = (max(Pr$FAge)/4)*3,
        y = 0.01,
        label = paste("k=", Coef[2]),
        hjust = "left"
      )
    
    name <- names(SSS[i])
    pf <- grid.arrange(p1, p2, nrow = 1, top = name)
    ggsave(
      plot = pf,
      path = nwpath,
      filename = paste(name, ".jpg"),
      units = "cm",
      width = 15,
      height = 10
    )
    
    
    #test.nlsResiduals(nlsResiduals(Exp1))
    
    
    mypath <-
      file.path(cidr, nwpath, paste(name, ".RES", ".jpg", sep = ""))
    
    jpeg(file = mypath)
    
    plot(nlsResiduals(Exp1))
    title(main = name)
    dev.off()
  }
  
  write.csv(Tfits,
            file.path(nwpath, "Tfits.csv"),
            sep = ";",
            dec = ".")
  return(Tfits)
} #function to estimate production decay models over oc acc mass per core. Needs df with
# columns: c("Core", "Ecosystem", "FAge", "Corg", "Corg.M") and a path to save outputs

OCModel(DataAM_lt, nwpath="Decay2023/Prueba")

#eliminate outlayers: Sm_004
#eliminate empty cores (no model): Sg_084,Sg_088;Sg_317, Sg_472, Sg_476, Sg_478, Sg_480, Sg_495, Sg_496, Sg_497 
#eliminate some cores after visual check, we eliminate: Mg_023, Sg_316, Sg_321, Sg_332, Sm_068, Sm_069, Sm_092, Sm_097, Sm_105
TfitsM_DEC <- TfitsM_DEC[-c(6,18, 19, 45, 46, 47, 49, 50:53, 55:57, 59, 72, 73, 78, 80, 81 ), ]



ggplot(TfitsM_DEC, aes(x = k)) +
  geom_histogram()

shapiro.test(TfitsM_DEC$k) #(>0.05 normal, <0.05 no normal)

# check if time frame has an effect 

ggplot(TfitsM_DEC, aes(Max.Age, k))+
  geom_point()

# results distribution

#get coordinates from B dataframe

CoordR<-subset(B, B$Core %in% TfitsM_DEC$ID==TRUE)

CoordR %>%
  ggplot() + ggtitle("Estimated k distribution") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
             pch = 21,
             size = 1.8) +
  coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
  scale_fill_manual(values = c( "blue", "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(path = Folder,
       filename =  "Estimated k (less_150) map.jpg",
       units = "cm",
       width = 20,
       height = 10
)

# only those between 100 and 150 time frame

TfitsM_DEC_100_150<-subset(TfitsM_DEC,TfitsM_DEC$Max.Age > 100 )

ggplot(TfitsM_DEC_100_150, aes(x = k)) +
  geom_histogram()

shapiro.test(TfitsM_DEC_100_150$k) #(>0.05 normal, <0.05 no normal)



CoordR<-subset(B, B$Core %in% TfitsM_DEC_100_150$ID==TRUE)

CoordR %>%
  ggplot() + ggtitle("Estimated k (100-150 years) distribution") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
             pch = 21,
             size = 2) +
  coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
  scale_fill_manual(values = c( "blue", "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(path = Folder,
       filename =  "Estimated k.Pb (100-150) map.jpg",
       units = "cm",
       width = 20,
       height = 10
)

# Comparison between dating methods  (same core, dif dating method)-------------------------------------

# First we estimate the correlation with depth of data.frame with those cores with both Pb and C. We use combined acc.mass-age model

ADT <- TPbC[, c("Core", "Ecosystem", "Age", "Corg")]

X <- split(ADT, ADT$Core)

TDT <- data.frame(
  ID = character(),
  Ecosystem = character(),
  C_rho = numeric(),
  C_p = numeric(),
  C_Gr = character()
)

for (i in 1:length(X)) {
  TDT[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(ADT)
  TDT[i, 2] <- Data[1, which(colnames(Data) == "Ecosystem")]
  #Data<-na.omit(Data)
  cor <- cor.test(x = Data[, 3], y = Data[, 4], method = "spearman")
  TDT[i, 3] <- cor$estimate
  TDT[i, 4] <- cor$p.value
  
}
##### from spearman correlations we discriminate 3 groups: NT, no trend with depth; DEC, decrease with depth; INC, increase with depth
#lower keys: IDs, capital keys: data frame with core data
# return data frames per group (NT,DEC and INC) and add trend tipe to DT


Nt <- TDT$ID[TDT$C_p > 0.5]
Trend <- TDT[!(TDT$C_p > 0.5), ]
Dec <- Trend$ID[Trend$C_rho < 0]
Inc <- Trend$ID[Trend$C_rho > 0]

TNT <- TPb[is.element(TPb$Core, Nt), ]
TDEC <- TPb[is.element(TPb$Core, Dec), ]
TINC <- TPb[is.element(TPb$Core, Inc), ]

for (i in 1:nrow(TDT)) {
  if (is.element(TDT[i, 1], Nt)) {
    TDT[i, 5] <- "NT"
  }
  else if (is.element(TDT[i, 1], Dec)) {
    TDT[i, 5] <- "DEC"
  }
  else {
    TDT[i, 5] <- "INC"
  }
}


# nlm 

Data <-
  as.data.frame(TDEC[, c("Core", "Ecosystem","FAge", "Corg", "Corg.M")])


Data$Corg.M <- as.numeric(Data$Corg.M)

SSS <- split(Data, Data$Core)


TfitsM_DEC <- data.frame(
  ID = character(),
  Ecosystem = character(),
  P = numeric(),
  k = numeric(),
  P.Pb = numeric(),
  k.Pb = numeric(),
  P.C = numeric(),
  k.C = numeric()
)

cidr <- getwd()
nwpath <- "Decay2023/AjustesComp"
dir.create(file.path(cidr, nwpath), recursive = TRUE)

for (i in 1:length(SSS)) {
  TfitsM_DEC[i, 1] <- names(SSS[i])
  Pr <- as.data.frame(SSS[i])
  colnames(Pr) <- list("Core", "Ecosystem", "Age", "Corg",  "Corg.M")
  TfitsM_DEC[i, 2] <- Pr[1, which(colnames(Data) == "Ecosystem")]
  
  skip_to_next <- FALSE
  
  # combined acc.mass-age model decay rates
  
  tryCatch(
    Exp1 <-
      nls(
        Corg.M ~ (p / k) * (1 - (exp(-k * Age))),
        data = Pr,
        start = list(p = 0.01, k = 0.03)
      )
    ,
    error = function(e) {
      skip_to_next <<- TRUE
    }
  )
  
  if (skip_to_next) {
    next
  }
  
  Func <-
    fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age))),
             data = Pr,
             start = list(p = 0.01, k = 0.03))
  
  Coef <- coef(Exp1)
  
  finales <- as.list(coef(Exp1))
  TfitsM_DEC[i, 3] <- as.numeric(finales[1])
  TfitsM_DEC[i, 4] <- as.numeric(finales[2])
  
  # Pb acc.mass-age model decay rates
  
  tryCatch(
    Exp1 <-
      nls(
        Corg.M ~ (p / k) * (1 - (exp(-k * Age.Pb))),
        data = Pr,
        start = list(p = 0.01, k = 0.03)
      )
    ,
    error = function(e) {
      skip_to_next <<- TRUE
    }
  )
  
  if (skip_to_next) {
    next
  }
  
  Func <-
    fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age.Pb))),
             data = Pr,
             start = list(p = 0.01, k = 0.03))
  
  Coef <- coef(Exp1)
  
  finales <- as.list(coef(Exp1))
  TfitsM_DEC[i, 5] <- as.numeric(finales[1])
  TfitsM_DEC[i, 6] <- as.numeric(finales[2])
  
  # 14C acc.mass-age model decay rates
  
  tryCatch(
    Exp1 <-
      nls(
        Corg.M ~ (p / k) * (1 - (exp(-k * Age.C))),
        data = Pr,
        start = list(p = 0.01, k = 0.03)
      )
    ,
    error = function(e) {
      skip_to_next <<- TRUE
    }
  )
  
  if (skip_to_next) {
    next
  }
  
  Func <-
    fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age.C))),
             data = Pr,
             start = list(p = 0.01, k = 0.03))
  
  Coef <- coef(Exp1)
  
  finales <- as.list(coef(Exp1))
  TfitsM_DEC[i, 7] <- as.numeric(finales[1])
  TfitsM_DEC[i, 8] <- as.numeric(finales[2])
  
}

ggplot(TfitsM_DEC, aes(k,k.Pb))+
  geom_point()+
  geom_point(aes(k,k.C), color="red")+
  xlim(0,0.05)+ylim(0,0.05)+
  geom_abline() 


#  Comparison between dating methods  (all cores) -------------------------

Data <-
  as.data.frame(TAll[, c("Core", "Ecosystem","Age", "Age.Pb", "Age.C", "Corg", "Corg.M")])


Data$Corg.M <- as.numeric(Data$Corg.M)

SAll <- split(Data, Data$Core)


TfitsM_All <- data.frame(
  ID = character(),
  Ecosystem = character(),
  P = numeric(),
  k = numeric(),
  P.Pb = numeric(),
  k.Pb = numeric(),
  P.C = numeric(),
  k.C = numeric()
)


for (i in 1:length(SAll)) {
  TfitsM_All[i, 1] <- names(SAll[i])
  Pr <- as.data.frame(SAll[i])
  colnames(Pr) <- list("Core", "Ecosystem","Age", "Age.Pb", "Age.C", "Corg", "Corg.M")
  TfitsM_All[i, 2] <- Pr[1, which(colnames(Pr) == "Ecosystem")]
  
  Pr1 <- Pr %>% filter(Age < 150)
  Pr2 <- Pr %>% filter(Age.Pb < 150)
  Pr3 <- Pr %>% filter(Age.C < 150)
  
  if ((count(!is.na(Pr1$Age)))>5) {
    
    cor <- cor.test(x = Pr1[,which(colnames(Pr1) == "Age")], y = Pr1[, which(colnames(Pr1) == "Corg")], method = "spearman")
    
    if (cor$p.value < 0.5 & cor$estimate < 0) {
      
      skip_to_next <- FALSE
      
      tryCatch(
        Exp1 <-
          nls(
            Corg.M ~ (p / k) * (1 - (exp(-k * Age))),
            data = Pr1,
            start = list(p = 0.01, k = 0.03)
          )
        ,
        error = function(e) {
          skip_to_next <<- TRUE
        }
      )
      
      if (skip_to_next) {
        next
      }
      
      Func <-
        fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age))),
                 data = Pr1,
                 start = list(p = 0.01, k = 0.03))
      
      Coef <- coef(Exp1)
      
      finales <- as.list(coef(Exp1))
      TfitsM_All[i, 3] <- as.numeric(finales[1])
      TfitsM_All[i, 4] <- as.numeric(finales[2])
      
      
      
    } }
  
  else { if (count(!is.na(Pr2$Age.Pb))>5) {
    cor <- cor.test(x = Pr2[,which(colnames(Pr2) == "Age.Pb")], y = Pr2[, which(colnames(Pr2) == "Corg")], method = "spearman")
      if (cor$p.value < 0.5 & cor$estimate < 0) {
        
        skip_to_next <- FALSE
        
        tryCatch(
          Exp1 <-
            nls(
              Corg.M ~ (p / k) * (1 - (exp(-k * Age.Pb))),
              data = Pr2,
              start = list(p = 0.01, k = 0.03)
            )
          ,
          error = function(e) {
            skip_to_next <<- TRUE
          }
        )
        
        if (skip_to_next) {
          next
        }
        
        Func <-
          fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age.Pb))),
                   data = Pr2,
                   start = list(p = 0.01, k = 0.03))
        
        Coef <- coef(Exp1)
        
        finales <- as.list(coef(Exp1))
        TfitsM_All[i, 5] <- as.numeric(finales[1])
        TfitsM_All[i, 6] <- as.numeric(finales[2])
        
      }}
    
    else { if (count(!is.na(Pr3$Age.C))>5) {
        cor <- cor.test(x = Pr3[,which(colnames(Pr3) == "Age.C")], y = Pr3[, which(colnames(Pr3) == "Corg")], method = "spearman")
        
        if (cor$p.value < 0.5 & cor$estimate < 0) {
          
          skip_to_next <- FALSE
          
          tryCatch(
            Exp1 <-
              nls(
                Corg.M ~ (p / k) * (1 - (exp(-k * Age.C))),
                data = Pr3,
                start = list(p = 0.01, k = 0.03)
              )
            ,
            error = function(e) {
              skip_to_next <<- TRUE
            }
          )
          
          if (skip_to_next) {
            next
          }
          
          Func <-
            fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age.C))),
                     data = Pr3,
                     start = list(p = 0.01, k = 0.03))
          
          Coef <- coef(Exp1)
          
          finales <- as.list(coef(Exp1))
          TfitsM_All[i, 7] <- as.numeric(finales[1])
          TfitsM_All[i, 8] <- as.numeric(finales[2])
        
      } }}}}


ggplot(TfitsM_All,aes(Ecosystem, k))+
  #geom_boxplot()+
  geom_jitter()+
  geom_jitter(aes(Ecosystem, k.Pb), color="red")+
  geom_jitter(aes(Ecosystem, k.C), color="blue")+
  ylim(-0.1,0.15)


#  Comparison between ecosystems  (all cores) -------------------------





#######################
### long term decay ###
#######################


ADT_lt <- TAllF[, c("Core", "Ecosystem","Min.Depth","Max.Depth", "FAge", "Corg")]

X <- split(ADT_lt, ADT_lt$Core)

TDT_lt <- data.frame(
  ID = character(),
  Ecosystem = character(),
  C_rho = numeric(),
  C_p = numeric(),
  C_Gr = character()
)



for (i in 1:length(X)) {
  TDT_lt [i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(ADT_lt)
  TDT_lt [i, 2] <- Data[1, which(colnames(Data) == "Ecosystem")]
  #Data<-na.omit(Data)
  if (nrow(Data)<3) next
  cor <- cor.test(x = Data$FAge, y = Data$Corg, method = "spearman")
  TDT_lt [i, 3] <- cor$estimate
  TDT_lt [i, 4] <- cor$p.value
  
}


##### from spearman correlations we discriminate 3 groups: NT, no trend with depth; DEC, decrease with depth; INC, increase with depth
#lower keys: IDs, capital keys: data frame with core data
# return data frames per group (NT,DEC and INC) and add trend tipe to DT


Nt <- TDT_lt $ID[TDT_lt$C_p > 0.5]
Trend <- TDT_lt[!(TDT_lt$C_p > 0.5), ]
Dec <- Trend$ID[Trend$C_rho < 0]
Inc <- Trend$ID[Trend$C_rho > 0]

TNT_lt <- TAllF[is.element(TAllF$Core, Nt), ]
TDEC_lt <- TAllF[is.element(TAllF$Core, Dec), ]
TINC_lt <- TAllF[is.element(TAllF$Core, Inc), ]

for (i in 1:nrow(TDT_lt)) {
  if (is.element(TDT_lt[i, 1], Nt)) {
    TDT_lt[i, 5] <- "NT"
  }
  else if (is.element(TDT_lt[i, 1], Dec)) {
    TDT_lt[i, 5] <- "DEC"
  }
  else {
    TDT_lt[i, 5] <- "INC"
  }
}



#### Count cores per group and get the percentages

library(dplyr)
NGr <- TDT_lt %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

#### plot per groups to visual check #####


ggplot(TNT_lt, aes(FAge, Corg)) + xlab("Time (years)") + ylab("Organic carbon (g cm-3)") +
  geom_point(aes(FAge, Corg)) +
  geom_line(aes(FAge, Corg)) +
  facet_wrap( ~ Core, ncol = 5 , scales = "free") +
  coord_flip() +
  scale_x_reverse() +
  theme_light()

ggsave(
  path = Folder,
  filename = 'TNT_lt.jpg',
  width = 20,
  height = 30,
  units = 'cm'
)

ggplot(TDEC_lt, aes(FAge, Corg)) + xlab("Time (years)") + ylab("Organic carbon (g cm-3)") +
  geom_point(aes(FAge, Corg)) +
  geom_line(aes(FAge, Corg)) +
  facet_wrap( ~ Core, ncol = 5, scales = "free") +
  coord_flip() +
  scale_x_reverse() +
  theme_light()

ggsave(
  path = Folder,
  filename = 'TDEC_lt.jpg',
  width = 20,
  height = 50,
  units = 'cm',
  limitsize = FALSE
)


ggplot(TINC_lt, aes(FAge, Corg)) + xlab("Time (years)") + ylab("Organic carbon (g cm-3)") +
  geom_point(aes(FAge, Corg)) +
  geom_line(aes(FAge, Corg)) +
  facet_wrap( ~ Core, ncol = 5 , scales = "free") +
  coord_flip() +
  scale_x_reverse() +
  theme_light()

ggsave(
  path = Folder,
  filename = 'TINC_lt.jpg',
  width = 20,
  height = 30,
  units = 'cm'
)


### Long term NLM time-Accumulated Mass ####

#Ajuste de cores que decrecen con la profundidad

DataA_lt <-
  as.data.frame(TDEC_lt[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])


#Estimate organic carbon accumulated mass

estimate_h <- function(df = NULL) {
  
  # create individual data frames per each core
  
  df$Core <- factor(df$Core, levels=unique(df$Core))
  X<-split(df, df$Core)
  
  
  columns<-c("EMin","EMax","h")
  Fdf2 = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(Fdf2) = columns
  
  
  for(i in 1:length(X)) {
    
    Data<-as.data.frame(X[i])
    colnames(Data)<-colnames(df)
    
    #check if there is spaces between samples (e.g, first sample ends at 5 cm and next starts at 7)
    space<- c()
    
    for (j in 1:(nrow(Data)-1)) {
      
      # if there are no spaces between samples min and maximun depth of samples remain the same
      if (Data[j,which(colnames(Data)=="Max.Depth")] == Data[j+1,which(colnames(Data)=="Min.Depth")]) {
        space[j]<-FALSE} else {space[j]<-TRUE}}
    
    if (any(space==TRUE)) {
      # if there are spaces between samples it estimate the medium point between the maximum depth of the sample and the minimum
      #depth of the next sample and divide that distance between both samples
      Data <- cbind(Data, EMin=NA, EMax=NA)
      Data[1,"EMin"]<-0
      Data[nrow(Data),"EMax"]<-Data[nrow(Data),"Max.Depth"]
      for (j in 1:(nrow(Data)-1)) {
        if(space[j]==TRUE) {
          Data[j,"EMax"]<-Data[j,"Max.Depth"]+((Data[j+1,"Min.Depth"]-Data[j,"Max.Depth"])/2)
          Data[j+1,"EMin"]<-Data[j,"Max.Depth"]+((Data[j+1,"Min.Depth"]-Data[j,"Max.Depth"])/2)} else {
            Data[j,"EMax"]<-Data[j,"Max.Depth"]
            Data[j+1,"EMin"]<-Data[j+1,"Min.Depth"]}}
      
    }  else{
      Data <- cbind(Data, EMin=NA, EMax=NA)
      Data$EMin<-Data$Min.Depth
      Data$EMax<-Data$Max.Depth
      
    }
    
    Data <- cbind(Data, h=NA)
    
    #estimation of the thickness of the sample (h) from the new minimun and max depth of the sample
    
    Data<- Data |> dplyr::mutate (h = EMax-EMin)
    
    temp<-cbind(Data$EMin, Data$EMax, Data$h)
    colnames(temp)<-colnames(Fdf2)
    Fdf2<-rbind(Fdf2, temp)
    
  }
  Fdf<-cbind(df, Fdf2)
  
  return(Fdf)
}

DataA_lt<-estimate_h(DataA_lt)


#estimate carbon density and acc mass per sample

#organic carbon mass per sample

DataA_lt<- DataA_lt %>% mutate (OCg = DBD*(Corg/100)*h)

#Acc organic matter

DataAM_lt<- DataA_lt[0,]
DataAM_lt[1,]=NA  # ad a temporary new row of NA values
DataAM_lt[,'Corg.M'] = NA # adding new column, called for example 'new_column'
DataAM_lt = DataAM_lt[0,]

X<- split(DataA_lt, DataA_lt$Core)

for (i in 1:length(X)) {
  Data <- as.data.frame(X[i])
  colnames(Data)<-colnames(DataA_lt)
  
  Data <- cbind(Data, Corg.M=NA)
  
  Data[1,"Corg.M"]<-Data[1,"OCg"]
  
  for (j in 2:nrow(Data)){
    Data[j,"Corg.M"]<-Data[j,"OCg"]+Data[j-1,"Corg.M"]
  }
  
  DataAM_lt<-rbind(DataAM_lt,Data)}


#estimate model time acc mat 

DataAM_lt<-as.data.frame(DataAM_lt[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
SSS <- split(DataAM_lt, DataAM_lt$Core)


TfitsM_DEC_lt <- data.frame(
  ID = character(),
  Ecosystem = character(),
  P = numeric(),
  k = numeric(),
  Res.Normality = numeric(),
  Res.Autocorrelation = numeric(),
  Max.Age = numeric()
)

cidr <- getwd()
nwpath <- "Decay2023/Ajustes_lt"
dir.create(file.path(cidr, nwpath), recursive = TRUE)

for (i in 1:length(SSS)) {
  TfitsM_DEC_lt[i, 1] <- names(SSS[i])
  Pr <- as.data.frame(SSS[i])
  colnames(Pr) <- list("Core", "Ecosystem", "FAge", "Corg", "Corg.M")
  TfitsM_DEC_lt[i, 2] <- Pr[1, which(colnames(Data) == "Ecosystem")]
  TfitsM_DEC_lt[i, 7] <- max(Pr$FAge)
  
  skip_to_next <- FALSE
  
  tryCatch(
    Exp1 <-
      nls(
        Corg.M ~ (p / k) * (1 - (exp(-k * FAge))),
        data = Pr,
        start = list(p = 0.01, k = 0.03)
      )
    ,
    error = function(e) {
      skip_to_next <<- TRUE
    }
  )
  
  if (skip_to_next) {
    next
  }
  
  Func <-
    fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * FAge))),
             data = Pr,
             start = list(p = 0.01, k = 0.03))
  
  Coef <- coef(Exp1)
  
  finales <- as.list(coef(Exp1))
  TfitsM_DEC_lt[i, 3] <- as.numeric(finales[1])
  TfitsM_DEC_lt[i, 4] <- as.numeric(finales[2])
  
  
  fitY1 <- as.data.frame(c(1:max(Pr$FAge)))
  fitY1['new_col'] <- NA
  fitY1[, 2] <- Func(c(1:max(Pr$FAge)), Coef[1], Coef[2])
  colnames(fitY1) <- list("FAge", "Predict")
  
  p1 <-
    ggplot(Pr, aes(FAge, Corg)) + xlab("Age (years)") + ylab("Corg (g cm-3)") +
    geom_point() +
    geom_line() +
    coord_flip() +
    scale_x_reverse()
  
  p2 <-
    ggplot(fitY1, aes(FAge, Predict)) + xlab("Age (years)") + ylab("Corg Aumulated mass (g cm-2)") +
    geom_line(color = "blue") +
    geom_point(data = Pr, aes(FAge, Corg.M)) +
    coord_flip() +
    scale_x_reverse() +
    annotate(
      "text",
      x = max(Pr$FAge)/4,
      y = 0.01,
      label = "M=(p/k)*(1-exp(-k*Age)",
      hjust = "left"
    ) +
    annotate(
      "text",
      x = max(Pr$FAge)/3,
      y = 0.01,
      label = paste("p=", Coef[1]),
      hjust = "left"
    ) +
    annotate(
      "text",
      x = max(Pr$FAge)/2,
      y = 0.01,
      label = paste("k=", Coef[2]),
      hjust = "left"
    )
  
  name <- names(SSS[i])
  pf <- grid.arrange(p1, p2, nrow = 1, top = name)
  ggsave(
    plot = pf,
    path = nwpath,
    filename = paste(name, ".jpg"),
    units = "cm",
    width = 15,
    height = 10
  )
  
  
  #test.nlsResiduals(nlsResiduals(Exp1))
  
  
  mypath <-
    file.path(cidr, nwpath, paste(name, ".RES", ".jpg", sep = ""))
  
  jpeg(file = mypath)
  
  plot(nlsResiduals(Exp1))
  title(main = name)
  dev.off()
}

write.csv(TfitsM_DEC_lt,
          file.path(nwpath, "TfitsM_DEC_lt.csv"),
          sep = ";",
          dec = ".")


#eliminate empty cores (no model): Sg_032, Sg_038, Sg_042, Sg_043, Sg_057, Sg_058, Sg_083, Sg_109, Sg_120, Sg_121, Sg_170, 
#Sg_191, Sg_193, Sg_241, Sg_147, Sg_315, Sg_317, Sg_483, Sg_484, Sg_491, Sg_492, Sg_498 
#eliminate some cores after visual check, we eliminate: Sg_019, Sg_041, Sg_088, Sg_097, Sg_310, Sg_311, Sg_312, Sg_314, Sg_316,
#Sg_321, Sg_476, Sg_485, Sg_490, Sg_495, Sg_497, Sm_010, Sm_022, Sm_068, Sm_069, Sm_092, Sm_097, Sm_105
TfitsM_DEC_lt <- TfitsM_DEC_lt[-c(15,12,16,17,18,19,21,22,24,28,26,27,33:34,36,41:44,51,52:57,58,61,68,69,70:73,76,78,79,80,81,94,95,100,102,103), ]



ggplot(TfitsM_DEC_lt, aes(x = k)) +
  geom_histogram()

shapiro.test(TfitsM_DEC$k) #(>0.05 normal, <0.05 no normal)

# check if time frame has an effect 

ggplot(TfitsM_DEC_lt, aes(Max.Age, k))+
  geom_point()

# results distribution

#get coordinates from B dataframe

CoordR<-subset(B, B$Core %in% TfitsM_DEC_lt$ID==TRUE)

CoordR %>%
  ggplot() + ggtitle("Estimated k distribution") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
             pch = 21,
             size = 1.8) +
  coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
  scale_fill_manual(values = c( "blue", "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(path = Folder,
       filename =  "Estimated k (all) map.jpg",
       units = "cm",
       width = 20,
       height = 10
)






































LTB <- B
LTB <- LTB[!(is.na(LTB$Age.Pb.C) | LTB$Age.Pb.C == ""),]


#we delete those cores with not good age-depth models
LTB <- LTB[!(is.na(LTB$Core) | LTB$Core == "Sg_097"),]
LTB <- LTB[!(is.na(LTB$Core) | LTB$Core == "Sg_112"),]
LTB <- LTB[!(is.na(LTB$Core) | LTB$Core == "Sg_121"),]

#spearman correlation between time and Corg content


ADT <- LTB[, c("Core", "Age.Pb.C", "Corg")]

X <- split(ADT, ADT$Core)

TDT <- data.frame(
  ID = character(),
  C_rho = numeric(),
  C_p = numeric(),
  C_Gr = character()
)



for (i in 1:length(X)) {
  TDT[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  #Data<-na.omit(Data)
  cor <- cor.test(x = Data[, 2], y = Data[, 3], method = "spearman")
  TDT[i, 2] <- cor$estimate
  TDT[i, 3] <- cor$p.value
  
}


##### from spearman correlations we discriminate 3 groups: NT, no trend with depth; DEC, decrease with depth; INC, increase with depth
#lower keys: IDs, capital keys: data frame with core data
# return data frames per group (NT,DEC and INC) and add trend tipe to DT


Nt <- TDT$ID[TDT$C_p > 0.5]
T <- TDT[!(TDT$C_p > 0.5), ]
Dec <- T$ID[T$C_rho < 0]
Inc <- T$ID[T$C_rho > 0]

TNT <- LTB[is.element(LTB$Core, Nt), ]
TDEC <- LTB[is.element(LTB$Core, Dec), ]
TINC <- LTB[is.element(LTB$Core, Inc), ]

for (i in 1:nrow(TDT)) {
  if (is.element(TDT[i, 1], Nt)) {
    TDT[i, 4] <- "NT"
  }
  else if (is.element(TDT[i, 1], Dec)) {
    TDT[i, 4] <- "DEC"
  }
  else {
    TDT[i, 4] <- "INC"
  }
}



#### Count cores per group and get the percentages

library(dplyr)
NGr <- TDT %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

######### long term nlm OC acc mass-age ##########

Data <- TDEC[, c("Core", "Age.Pb.C", "Corg", "Corg.M")]

Data$M <- as.numeric(Data$Corg.M)

SSS <- split(Data, Data$Core)


TfitsM_DEC_long <- data.frame(
  ID = character(),
  Tframe = numeric(),
  P = numeric(),
  k = numeric()
)

cidr <- getwd()
nwpath <- "SG/ajustesDECLong"
dir.create(file.path(cidr, nwpath), recursive = TRUE)



for (i in 1:length(SSS)) {
  TfitsM_DEC_long[i, 1] <- names(SSS[i])
  Pr <- as.data.frame(SSS[i])
  colnames(Pr) <- list("Core", "Age", "Corg", "Corg.M")
  TfitsM_DEC_long[i, 2] <- max(Pr$Age)
  
  Exp1 <-
    nls(Corg.M ~ (p / k) * (1 - (exp(-k * Age))),
        data = Pr,
        start = list(p = 0.01, k = 0.03))
  Func <-
    fitModel(Corg.M ~ (p / k) * (1 - (exp(-k * Age))),
             data = Pr,
             start = list(p = 0.01, k = 0.03))
  
  Coef <- coef(Exp1)
  
  finales <- as.list(coef(Exp1))
  TfitsM_DEC_long[i, 3] <- as.numeric(finales[1])
  TfitsM_DEC_long[i, 4] <- as.numeric(finales[2])
  
  
  fitY1 <- as.data.frame(c(1:max(Pr$Age)))
  fitY1['new_col'] <- NA
  fitY1[, 2] <- Func(c(1:max(Pr$Age)), Coef[1], Coef[2])
  colnames(fitY1) <- list("Age", "Predict")
  
  p1 <-
    ggplot(Pr, aes(Age, Corg)) + xlab("Age (years)") + ylab("Corg (g cm-3)") +
    geom_point() +
    geom_line() +
    coord_flip() +
    scale_x_reverse()
  
  p2 <-
    ggplot(fitY1, aes(Age, Predict)) + xlab("Age (years)") + ylab("Corg Aumulated mass (g cm-3)") +
    geom_line(color = "blue") +
    geom_point(data = Pr, aes(Age, Corg.M)) +
    coord_flip() +
    scale_x_reverse() +
    annotate(
      "text",
      x = 60,
      y = 0.01,
      label = "M=(p/k)*(1-exp(-k*Age)",
      hjust = "left"
    ) +
    annotate(
      "text",
      x = 100,
      y = 0.01,
      label = paste("p=", Coef[1]),
      hjust = "left"
    ) +
    annotate(
      "text",
      x = 150,
      y = 0.01,
      label = paste("k=", Coef[2]),
      hjust = "left"
    )
  
  name <- names(SSS[i])
  pf <- grid.arrange(p1, p2, nrow = 1, top = name)
  ggsave(
    plot = pf,
    path = nwpath,
    filename = paste(name, ".jpg"),
    units = "cm",
    width = 15,
    height = 10
  )
  
  
  #test.nlsResiduals(nlsResiduals(Exp1))
  
  
  mypath <-
    file.path(cidr, nwpath, paste(name, ".RES", ".jpg", sep = ""))
  
  jpeg(file = mypath)
  
  plot(nlsResiduals(Exp1))
  title(main = name)
  dev.off()
  
  
}


write.csv(
  TfitsM_DEC_long,
  file.path(nwpath, "TfitsM_DEC_long.csv"),
  sep = ";",
  dec = "."
)

ggplot(TfitsM_DEC_long, aes(x = k)) +
  geom_histogram()

#######################
### summary results ###
#######################
## ajustes tras inspeccionar visualemnte....


File <- "Decay150.csv"

R <- read.csv(File,
              header = T,
              sep = ";",
              dec = ".")
R <- as.data.frame(R)

# map


WM <- map_data("world")


R %>%
  ggplot() + ggtitle("Decay rates sampling site") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  geom_point(aes(x = long, y = lat,  fill = Ecosystem),
             pch = 21,
             size = 2) +
  coord_sf(ylim = c(-40, 50), xlim = c(-90, 140), ) +
  scale_fill_manual(values = c("yellow", "blue", "orange", "green")) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave("Decay_map.jpg",
       units = "cm",
       width = 20,
       height = 10)

#summary

summary(R$k)

ggplot(R, aes(x = k)) +
  geom_histogram()


ggplot(R, aes(Tipo, k)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir = 'center',
               dotsize = 0.5)

shapiro.test(R$k) #(>0.05 normal, <0.05 no normal)

## Student's t-test  if normally distributed, wilcox if not

pairwise.wilcox.test(R$k, R$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


library(reshape2)

R2 <- dcast(R, k ~ Ecosystem)

t.test(R[c(19:30), 3], R[c(31:34), 3])

### plots ###

File <- "DecayS.csv"

P <- read.csv(File,
              header = T,
              sep = ";",
              dec = ".")
P <- as.data.frame(P)

P1 <- P[c(1:8), ]
P2 <- P[-c(1:8), ]


ggplot() + ggtitle("Decay rate by time frame in Posidonia spp. meadows") + xlab("Time frame (years)") + ylab("Decay rate (yr-1)") +
  geom_point(aes(P1$Tframe, P1$k), size = 3, shape = 16) +
  geom_point(aes(P2$Tframe, P2$k), size = 4, shape = 17)

plot(P$Tframe, P$k)

### exponential model to predict k in Posidonia meadows


kchange <- function(Tframe, A, C)
  (A * exp(C * Tframe))

model1 <-
  nls(
    k ~ kchange(Tframe, myA, myC),
    data = P,
    start = list(myA = 0.04, myC = -0.00201)
  )


summary(model1)
plot(nlsResiduals(model1))

fitY1 <- kchange(c(1:2100), 0.0401436, -0.0027)

P1 <- P[c(1:5, 7:8), ]
P2 <- P[-c(1:8), ]

plot(P$Tframe,
     P$k,
     main = expression(paste(
       italic("Posidonia spp."),
       " organic carbon decay rate by time frame"
     )),
     xlab = "Time frame (yr)",
     ylab = "Decay rate (yr-1)")
points(P1$Tframe, P1$k, pch = 16)
points(P2$Tframe, P2$k)
lines(c(1:2100), fitY1, col = "blue")
text(1000, 0.025, cex = 1.5, expression(y == 0.0401436 * e ** (-0.0027501 *
                                                                 x)))



0.0401436 * exp(-0.0027501 * 100)
