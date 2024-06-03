################################################################################
# Article: Modelling decay rates of soil organic carbon in blue carbon ecosystems 
# 1/1
# Author: Nerea Piñeiro-Juncal (https://github.com/NPJuncal)
# V 3.3
# Year: 2024
################################################################################

setwd("C:/Users/npjun/Dropbox/Seagrasses/Degradacion anaerobia_OSCAR/SOC_Decay")

library(ggplot2)
library(maps)
library(tidyr)
library(magrittr)
library(mosaic)
library(gridExtra)
library(nlstools)
library(dplyr)
library(MASS) 
library(reshape2) 
library(reshape) 
library(psych)



# Loading the dataset ------------------------------------------------------

File <- "Data/Cores.csv"

Cores <- read.csv(File,
                  header = T,
                  sep = ";",
                  dec = ".")
Cores <- as.data.frame(Cores)

# check numeric variables in dataset

is.numeric(Cores$DBD)
is.numeric(Cores$Corg)
is.numeric(Cores$Corg_.gcm3)
is.numeric(Cores$d13C)
is.numeric(Cores$Mud)

Cores$Max.Depth<-as.numeric(Cores$Max.Depth)

#clean the rows without DBD and OC data
Cores$DBD <- as.numeric(Cores$DBD)
Cores <- Cores[!is.na(Cores$DBD),]
Cores$Corg <- as.numeric(Cores$Corg)
Cores <- Cores[!is.na(Cores$Corg),]

# estimation of carbon density for those samples where there are bulk density and oc percentage available

Cores$Corg_.gcm3 <- as.numeric(Cores$Corg_.gcm3)

for (i in 1:nrow(Cores)) {
  
  if (is.na(Cores[i, "Corg_.gcm3"]) & !is.na(Cores[i, "DBD"]) & !is.na(Cores[i, "Corg"])) {
    
    Cores[i, "Corg_.gcm3"] <- Cores[i, "DBD"] * (Cores[i, "Corg"]/100)}}




### we create a folder to save the results

Folder = "Decay2023"
dir.create(Folder)


### We select only vegetated or un-vegetated (Bare) cores

unique(Cores[, 5])

B <- subset(Cores, Cores$V.vs.B == "Vegetated" | Cores$V.vs.B == "vegetated")
unique(B[, 5])


B$Corg <- as.numeric(B$Corg)
B$Mud <- as.numeric(B$Mud)
B$DBD <- as.numeric(B$DBD)


# we only keep cores with more than 4 samples


temp <- B[0,]


X<-split(B, B$Core)

for (i in 1:length(X)) {
  
  data<-X[[i]]
  
  if(nrow(data)>4) {
    
    temp<-rbind(temp, data)}}

length(unique(temp$Core))

B<-temp

#


length(unique(B$Core))
SingleCore<-B[!duplicated(B$Core),]

table(SingleCore$Ecosystem)

write.csv(SingleCore[,c(1,3,9,10)],
          file.path(Folder, "singlecore.csv"),
          sep = ";",
          dec = ".")


# Sampling sites Map (Figure 4) -------------------------------------------



#load a world map
WM <- map_data("world")

global<- B %>%
  ggplot() + ggtitle("Sampling sites") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
  coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
  scale_fill_manual(values = c("blue",  "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5))


nam<-B %>%
  ggplot() + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
  coord_sf(xlim = c(-150, -50), ylim = c(-20, 80)) +
  scale_fill_manual(values = c("blue",  "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")


eu<-B %>%
  ggplot()  + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
  coord_sf(xlim = c(-10, 50), ylim = c(20, 60)) +
  scale_fill_manual(values = c("blue",  "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")


aus<-B %>%
  ggplot()  + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = long, y = lat))+
  geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
  coord_sf(xlim = c(110, 155), ylim = c(-40, -5))  +
  scale_fill_manual(values = c("blue",  "green","orange")) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")

ssp<-grid.arrange(global, nam, eu, aus, 
                  layout_matrix = rbind(c(1, 1, 1),
                                        c(2, 3, 4)))


ggsave(
  plot = ssp,
  path = Folder,
  filename =  "Sampling sites.jpg",
  units = "cm",
  width = 20,
  height = 15
)

# Average and median OC content per core (full length and top 25 c --------

ADT <- B[, c("Core", "Ecosystem", "Max.Depth", "Corg", "Mud", "d13C")]

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
  Av_13C = numeric(),
  SE_13C = numeric(),
  M_13C = numeric(),
  Av_C_25 = numeric(),
  SE_C_25 = numeric(),
  M_C_25 = numeric(),
  Av_Mud_25 = numeric(),
  SE_Mud_25 = numeric(),
  M_Mud_25 = numeric(),
  Av_13C_25 = numeric(),
  SE_13C_25 = numeric(),
  M_13C_25 = numeric()
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
  CM[i, 9] <- mean(Data[, which(colnames(Data) == "d13C")], na.rm=TRUE)
  CM[i, 10] <- std(Data[, which(colnames(Data) == "d13C")])
  CM[i, 11] <- median(Data[, which(colnames(Data) == "d13C")], na.rm=TRUE)
  
  Data25 <-
    Data %>% filter(Data[, which(colnames(Data) == "Max.Depth")] <= 25)
  CM[i, 12] <- mean(Data25[, which(colnames(Data25) == "Corg")], na.rm=TRUE)
  CM[i, 13] <- std(Data25[, which(colnames(Data25) == "Corg")])
  CM[i, 14] <- median(Data25[, which(colnames(Data) == "Corg")], na.rm=TRUE)
  CM[i, 15] <- mean(Data25[, which(colnames(Data25) == "Mud")], na.rm=TRUE)
  CM[i, 16] <- std(Data25[, which(colnames(Data25) == "Mud")])
  CM[i, 17] <- median(Data[, which(colnames(Data) == "Mud")], na.rm=TRUE)
  CM[i, 18] <- mean(Data25[, which(colnames(Data25) == "d13C")], na.rm=TRUE)
  CM[i, 19] <- std(Data25[, which(colnames(Data25) == "d13C")])
  CM[i, 20] <- median(Data[, which(colnames(Data) == "d13C")], na.rm=TRUE)
}


write.csv(CM,
          file.path(Folder, "AvMdC.csv"),
          sep = ";",
          dec = ".")


# Sediment accretion rate -------------------------------------------------

File <- "Data/Acc_Mass-Age_F.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ";",
                  dec = ".")
Dates <- as.data.frame(Dates)

# estimates SAR last 150 years 

SAR <- data.frame(
  ID = character(),
  Ecosystem = character(),
  SAR_Age = numeric(),
  SAR_Pb = numeric(),
  SAR_C = numeric(),
  SAR_25 = numeric())

X<-split(Dates, Dates$Core)

for (i in 1:length(X)) {
  SAR[i, 1] <- names(X[i])
  Data <- as.data.frame(X[i])
  colnames(Data) <- colnames(Dates)
  Eco <- substr(Data[1,"Core"], 1, 2)
  
  if (Eco == "Mg") { SAR[i, 2] <- "Mangrove"}
  if (Eco == "Sg") { SAR[i, 2] <- "Seagrass"}
  if (Eco == "Sm") { SAR[i, 2] <- "Tidal Marsh"}
  
  
  #correlation depth-age to predict depth at 150 yr old
  DataAge <- Data[!is.na(Data$Age),]
  DataAge<- DataAge[c(1:(length(which(DataAge$Age <=150)))),]
  
  if (nrow(DataAge)>2) {SAR[i, 3] <- max(DataAge$Depth)/max(DataAge$Age)}
  
  DataPb <- Data[!is.na(Data$Age.Pb),]
  DataPb<- DataPb[c(1:(length(which(DataPb$Age.Pb <=150)))),]
  
  if (nrow(DataPb)>2) {SAR[i, 4] <- max(DataPb$Depth)/max(DataPb$Age.Pb)}
  
  DataC <- Data[!is.na(Data$Age.C),]
  DataC<- DataC[c(1:(length(which(DataC$Age.C <=150)))),]
  
  if (nrow(DataC)>2) {SAR[i, 5] <- max(DataC$Depth)/max(DataC$Age.C)}
  
  if (nrow(DataAge)>2) {
    Data25 <- Data[!is.na(Data$Age),]
    Data25<- Data25[c(1:(length(which(Data25$Depth <=25)))),]
    
    if (nrow(Data25)>2) {SAR[i, 6] <- max(Data25$Depth)/max(Data25$Age)}
  
    } else {
    
    Data25 <- Data[!is.na(Data$Age.Pb),]
    Data25<- Data25[c(1:(length(which(Data25$Depth <=25)))),]
    
    if (nrow(Data25)>2) {SAR[i, 6] <- max(Data25$Depth)/max(Data25$Age.Pb)}}
}

SAR<-SAR[-1,]
SAR$SAR<-"NA"
SAR$SAR<-as.numeric(SAR$SAR)

for (i in 1:nrow(SAR)) {
  
  if (is.na(SAR[i,"SAR_Age"]) == FALSE) {SAR[i,"SAR"]<-SAR[i,"SAR_Age"]} 
  
  else {
    
    if (is.na(SAR[i,"SAR_Pb"]) == FALSE) {SAR[i,"SAR"]<-SAR[i,"SAR_Pb"]}
    
    else {if (is.na(SAR[i,"SAR_C"]) == FALSE) {SAR[i,"SAR"]<-SAR[i,"SAR_C"]}}}}


ggplot(SAR, aes(Ecosystem, SAR)) +
  geom_boxplot() +
  geom_jitter()

pairwise.wilcox.test(SAR$SAR, SAR$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

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

#### Count cores per group and get the percentages

NGr <- DT %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

DT %>% group_by(Ecosystem, C_Gr) %>% count()
NGrSg <- subset(DT, Ecosystem == "Seagrass") %>% group_by(C_Gr) %>% count()
NGrSg %>% mutate(proc = ((n * 100) / sum(NGrSg[, 2])))

NGrSm <-
  subset(DT, Ecosystem == "Tidal Marsh") %>% group_by(C_Gr) %>% count()
NGrSm %>% mutate(proc = ((n * 100) / sum(NGrSm[, 2])))

NGrMg <- subset(DT, Ecosystem == "Mangrove") %>% group_by(C_Gr) %>% count()
NGrMg %>% mutate(proc = ((n * 100) / sum(NGrMg[, 2])))


# Figures and tables tendency with depth ----------------------------------


### Figure 1 ###



DT2<-DT
DT2$C_Gr <- recode(DT2$C_Gr, DEC = 'Decreasing',
                   INC  = 'No decreasing',
                   NT = 'No decreasing')



DT2<-cbind(CM, DT2)
DT2<-DT2[,c(1:20, 25)]
DT2 <-merge(DT2, SAR[,c(1,3:7)], by = "ID", all = TRUE)

cor.test(DT2$Av_Mud_25, DT2$Av_13C_25)
plot(DT2$Av_Mud_25, DT2$Av_13C_25)

### data from BIO-Oracle ##

File <- "Data/bio_oracle.csv"

oracle <- read.csv(File,
                   header = T,
                   sep = ";",
                   dec = ".")
oracle <- as.data.frame(oracle)


DT3 <-merge(DT2, oracle[c(1,7:12)], by = "ID", all = TRUE)
DT3<-subset(DT3, !is.na(DT3$C_Gr))

# seagrass specie and trend -----------------------------------------------

  Sg_DT<-subset(DT, Ecosystem=="Seagrass")
  Sg_DT$Life.form<-"NA"
  
  for (i in 1:nrow(Sg_DT)) {
    
    life_form<- unique(SingleCore[c(which(SingleCore$Core==Sg_DT[i,which(colnames(Sg_DT)=="ID")])),which(colnames(SingleCore)=="Life.form")])
    Sg_DT[i,which(colnames(Sg_DT)=="Life.form")]<- life_form
    
  }
  
  
  # percentage of dec in persistent and another life forms
  
  #persistent
  temp <- subset(Sg_DT, Life.form=="Persistent") %>% group_by(C_Gr) %>% count()
  temp %>% mutate(proc = ((n * 100) / sum(temp[, 2])))

  #no persistent
  temp <- subset(Sg_DT, !Life.form=="Persistent") %>% group_by(C_Gr) %>% count()
  temp %>% mutate(proc = ((n * 100) / sum(temp[, 2])))

  # percentage of life.forms in dec trends
  
  temp <- subset(Sg_DT, C_Gr=="DEC") %>% group_by(Life.form) %>% count()
  temp %>% mutate(proc = ((n * 100) / sum(temp[, 2])))
  
  # percentage of life.forms in dec trends
  
  temp <- subset(Sg_DT, !C_Gr=="DEC") %>% group_by(Life.form) %>% count()
  temp %>% mutate(proc = ((n * 100) / sum(temp[, 2])))
  
  
  
  ##
  
  #Seagrass meadows
  DT2Sg<-subset(DT2, Ecosystem=="Seagrass")
  
  pairwise.wilcox.test(DT2Sg$Av_13C, DT2Sg$Life.form,
                       p.adjust.method = "BH") # are significantly different (p < 0.05)
  
  

# Figure 1 ----------------------------------------------------------------


#Seagrass meadows
DT2Sg<-subset(DT2, Ecosystem=="Seagrass")



cor.test(DT2Sg$Av_C_25, DT2Sg$Av_13C_25, method="spearman")
plot(DT2Sg$Av_C_25, DT2Sg$Av_13C_25)

#organic carbon  

SC<-ggplot(DT2Sg, aes(C_Gr, Av_C_25)) + ggtitle("Seagrass")+
  geom_boxplot()+
  geom_jitter(color="green4", alpha = 0.1)+
  scale_y_log10()+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = 0.003,
           label = table(subset(DT2Sg, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())





#mud  
SM<-ggplot(DT2Sg, aes(C_Gr, Av_Mud_25)) + 
  geom_boxplot()+
  geom_jitter(color="green4", alpha = 0.5)+
  ylim(-10,110)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = -7,
           label = table(subset(DT2Sg, !is.na(Av_Mud_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())



#SAR

SS<-ggplot(DT2Sg, aes(C_Gr, SAR_25)) + 
  geom_boxplot()+
  geom_jitter(color="green4", alpha = 0.5)+
  ylim(-0.15,1.1)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = -0.1,
           label = table(subset(DT2Sg, !is.na(SAR_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())



#d13C

d13S<-ggplot(DT2Sg, aes(C_Gr, Av_13C_25)) + ylab(expression(delta~"13C (‰)")) + ggtitle("Seagrass")+
  geom_boxplot()+
  geom_jitter(color="green4", alpha = 0.3)+
  ylim(-30,-5)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = -30,
           label = table(subset(DT2Sg, !is.na(Av_13C_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),)

 
#temperature  

ST<-ggplot(subset(DT3, Ecosystem=="Seagrass"), aes(C_Gr, Temperature)) +
  geom_boxplot()+
  geom_jitter(color="green4", alpha = 0.1)+
  ylim(2,35)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = 2.5,
           label = table(subset(DT2Sg, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())


#current velocity  

SCV<-ggplot(subset(DT3, Ecosystem=="Seagrass"), aes(C_Gr, current.velocity)) + 
  ylab(expression(paste("Current v. (m ",s^-1,")"))) + ggtitle("Seagrass")+
  geom_boxplot()+
  geom_jitter(color="green4", alpha = 0.3)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = -0.05,
           label = table(subset(DT2Sg, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  theme(plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank())





#Tidal marshes  
DT2Sm<-subset(DT2, Ecosystem=="Tidal Marsh")

#organic carbon  
TMC<-ggplot(DT2Sm, aes(C_Gr, Av_C_25)) + ggtitle("Tidal Marsh")+
  geom_boxplot()+
  geom_jitter(color="orange", alpha = 0.1)+
  scale_y_log10()+
  annotate("text",
           x = 1:length(table(DT2Sm$C_Gr)),
           y = 0.003,
           label = table(subset(DT2Sm, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())


#mud  
TMM<-ggplot(DT2Sm, aes(C_Gr, Av_Mud_25)) +
  geom_boxplot()+
  geom_jitter(color="orange", alpha = 0.5)+
  ylim(-10,110)+
  annotate("text",
           x = 1:length(table(DT2Sm$C_Gr)),
           y = -7,
           label = table(subset(DT2Sm, !is.na(Av_Mud_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())



#SAR 
TMS<-ggplot(DT2Sm, aes(C_Gr, SAR_25)) +
  geom_boxplot()+
  geom_jitter(color="orange", alpha = 0.5)+
  ylim(-0.15,1.1)+
  annotate("text",
           x = 1:length(table(DT2Sm$C_Gr)),
           y = -0.1,
           label = table(subset(DT2Sm, !is.na(SAR_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())



#d13C

d13T<-ggplot(DT2Sm, aes(C_Gr, Av_13C_25)) + 
  geom_boxplot()+
  geom_jitter(color="orange", alpha = 0.5)+
  annotate("text",
           x = 1:length(table(DT2Sm$C_Gr)),
           y = -30,
           label = table(subset(DT2Sm, !is.na(Av_13C_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())


#temperature  

TT<-ggplot(subset(DT3, Ecosystem=="Tidal Marsh"), aes(C_Gr, Temperature)) +
  geom_boxplot()+
  geom_jitter(color="orange", alpha = 0.1)+
  ylim(2,35)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = 2.5,
           label = table(subset(DT2Sg, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())



#Mangroves        
DT2Mg<-subset(DT2, Ecosystem=="Mangrove")
#organic carbon  
MGC<-ggplot(DT2Mg, aes(C_Gr, Av_C_25)) + ylab("log(OC%)") + ggtitle("Mangrove")+
  geom_boxplot()+
  geom_jitter(color="blue", alpha = 0.1)+
  annotate("text",
           x = 1:length(table(DT2Mg$C_Gr)),
           y = 0.003,
           label = table(subset(DT2Mg, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  scale_y_log10()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())


#mud  
TGM<-ggplot(DT2Mg, aes(C_Gr, Av_Mud_25)) + ylab("Mud%") +
  geom_boxplot()+
  geom_jitter(color="blue", alpha = 0.5)+
  ylim(-10,110)+
  annotate("text",
           x = 1:length(table(DT2Mg$C_Gr)),
           y = -7,
           label = table(subset(DT2Mg, !is.na(Av_Mud_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank())



#SAR 
TGS<-ggplot(DT2Mg, aes(C_Gr, SAR_25)) + ylab(expression(paste("Sed. acc. (cm ",yr^-1,")"))) +
  geom_boxplot()+
  geom_jitter(color="blue", alpha = 0.5)+
  ylim(-0.15,1.1)+
  annotate("text",
           x = 1:length(table(DT2Mg$C_Gr)),
           y = -0.1,
           label = table(subset(DT2Mg, !is.na(SAR_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank())



#d13C

d13M<-ggplot(DT2Mg, aes(C_Gr, Av_13C_25)) + ylab(expression(delta~"13C (‰)")) +
  geom_boxplot()+
  geom_jitter(color="blue", alpha = 0.5)+
  ylim(-30,-5)+
  annotate("text",
           x = 1:length(table(DT2Mg$C_Gr)),
           y = -30,
           label = table(subset(DT2Mg, !is.na(Av_13C_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank())


#temperature  

MT<-ggplot(subset(DT3, Ecosystem=="Mangrove"), aes(C_Gr, Temperature)) +  ylab("Sea Water Tª (\u00B0C)") +
  geom_boxplot()+
  geom_jitter(color="blue", alpha = 0.1)+
  ylim(2,35)+
  annotate("text",
           x = 1:length(table(DT2Sg$C_Gr)),
           y = 2.5,
           label = table(subset(DT2Sg, !is.na(Av_C_25))[,"C_Gr"]),
           col = "black")+
  theme(axis.title.x = element_blank())




todos_CM<-ggpubr::ggarrange(MGC, SC, TMC, TGM, SM, TMM, TGS, SS, TMS, MT, ST, TT, d13S, SCV,
                            labels = c("A","B", "C","D" , "E","F", "G","H", "I", "J", "K", "L", "M", "N"),
                            ncol=3, nrow= 5)


ggsave( plot = todos_CM,
        path = Folder,
        filename =  "C_Mud_Gr.jpg",
        units = "cm",
        width = 20,
        height = 25)    

    ## diff sig 
    
    #Seagrass
    
    pairwise.wilcox.test(DT2Sg$Av_C_25, DT2Sg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    pairwise.wilcox.test(DT2Sg$Av_Mud_25, DT2Sg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)  
    
    pairwise.wilcox.test(DT2Sg$SAR_25, DT2Sg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    pairwise.wilcox.test(subset(DT3, Ecosystem=="Seagrass")$Temperature, subset(DT3, Ecosystem=="Seagrass")$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05) 
    
    pairwise.wilcox.test(DT2Sg$Av_13C_25, DT2Sg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)  
    
    pairwise.wilcox.test(subset(DT3, Ecosystem=="Seagrass")$current.velocity, subset(DT3, Ecosystem=="Seagrass")$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)  
    
    
    #Tidal marshes  
    
    
    pairwise.wilcox.test(DT2Sm$Av_C_25, DT2Sm$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    pairwise.wilcox.test(DT2Sm$Av_Mud_25, DT2Sm$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)  
    
    pairwise.wilcox.test(DT2Sm$SAR_25, DT2Sm$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05) 

    
    pairwise.wilcox.test(subset(DT3, Ecosystem=="Tidal Marsh")$Temperature, subset(DT3, Ecosystem=="Tidal Marsh")$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)  
    
    
    #Mangroves        
    
    
    pairwise.wilcox.test(DT2Mg$Av_C_25, DT2Mg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    pairwise.wilcox.test(DT2Mg$Av_Mud_25, DT2Mg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    pairwise.wilcox.test(DT2Mg$SAR_25, DT2Mg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05) 

    
    pairwise.wilcox.test(subset(DT3, Ecosystem=="Mangrove")$Temperature, subset(DT3, Ecosystem=="Mangrove")$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)  
    
    
    
    # average per ecosystem
    
    std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
    
    aggregate(DT2[,c(12, 15, 18, 25)], list(DT2$Ecosystem), FUN=mean, na.rm=T)
    aggregate(DT2[,c(12, 15, 18, 25)], list(DT2$Ecosystem), FUN=std) 
    
    
    aggregate(DT3[,c(27:32)], list(DT3$Ecosystem), FUN=mean, na.rm=T)
    aggregate(DT3[,c(6:11)], list(DT3$Ecosystem), FUN=std) 
    
# correlation with time ---------------------------------------------------



### chronological models estimation

File <- "Data/Acc_Mass-Age_F.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ";",
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
    Data <- cbind(Data, temp[, c(3:6)])
    
    #check if there is at least 5 samples with dates
    if (count(!is.na(Data$Age))>5 | count(!is.na(Data$Age.Pb))>5 | count(!is.na(Data$Age.C))>5) {
    
    C <- rbind(C, Data)}
    
  } else
    next
}

length(unique(C$Core))

single_core_dates<-C[!duplicated(C$Core),]


table(single_core_dates$Ecosystem)

table(subset(single_core_dates,Ecosystem=="Seagrass")[,"Specie"])


TAll = filter(C, !is.na(Age) | !is.na(Age.Pb) | !is.na(Age.C))# df with core with any model
  
  length(unique(TAll$Core))
  # number of cores with age-depth models per ecosystem
  table((TAll[!duplicated(TAll$Core),])[,"Ecosystem"]) 

TPb = subset(C, !is.na(Age) | !is.na(Age.Pb) )# df with core with Pb model

  length(unique(TPb$Core))

TC = filter(C, !is.na(Age.C))# df with core with C model
  
  length(unique(TC$Core))


TPbandC = filter(C, !is.na(Age) | !is.na(Age.Pb))# df with core with C and Pb model
TPbandC = filter(TPbandC, !is.na(Age.C))
length(unique(TPbandC$Core))


#### Homogenize Age


# Decay models ---------------------------------------------------



# function to generate different df depending on their tendency with time (DEC, INC and NT), generate figures per group 
# and return a dfs: TDT, TDEC, TINC and TNT. Needs a df with columns: c("Core", "Ecosystem","Min.Depth","Max.Depth", "FAge", "Corg"
# a max age (MA) and a name for the plots 

tendency<- function (df, pnames) {
  
  ADT <- df[, c("Core", "Ecosystem","Min.Depth","Max.Depth", "FAge", "Corg")]
  
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
  
  
  df_list<-list(TDT, TDEC, TINC, TNT)
  names(df_list)<-c("Sp_TDT", "TDEC", "TINC", "TNT")
  return(df_list)
}


estimate_h <- function(df = NULL,
                       core = "core",
                       mind = "mind",
                       maxd = "maxd") {
  
  # class of the dataframe or tibble
  if (!inherits(df, "data.frame")) {
    stop("The data provided must be a tibble or data.frame")
  }
  
  # name of the columns
  if (!core %in% colnames(df)) {stop("There must be a variable with 'core'")}
  if (!mind %in% colnames(df)) {stop("There must be a variable with 'mind'")}
  if (!maxd %in% colnames(df)) {stop("There must be a variable with 'maxd'")}
  
  # class of the columns
  if (!is.numeric(df[[mind]])) {stop("'mind' data must be class numeric")}
  if (!is.numeric(df[[maxd]])) {stop("'maxd' data must be class numeric")}
  
  #check for NAs in depth columns
  if (sum(is.na(df[[mind]])) > 0) {stop("Samples minimun depth column has NAs, please check")}
  if (sum(is.na(df[[maxd]])) > 0) {stop("Samples maximun depth column has NAs, please check")}
  
  # create variables with working names with the data in the columns specified by the user
  df_r <- df
  df_r$core_r <- df_r[[core]]
  df_r$mind_r <- df_r[[mind]]
  df_r$maxd_r <- df_r[[maxd]]
  
  # create individual data frames per each core
  df_r$core_r <- factor(df_r$core_r, levels = unique(df_r$core_r))
  x <- split(df_r, df_r$core_r)
  
  estimate_depth <- function(df, j) {
    df[j + 1, "emin"] <- df[j, "maxd_r"] + ((df[j + 1, "mind_r"] - df[j, "maxd_r"]) / 2)
    df[j, "emax"] <- df[j, "maxd_r"] + ((df[j + 1, "mind_r"] - df[j, "maxd_r"]) / 2)
    df[1, "emin"] <- 0
    df[nrow(df), "emax"] <- df[nrow(df), "maxd_r"]
    return(df)
  }
  
  estimate_height <- function(df) {
    data <- as.data.frame(df)
    colnames(data) <- colnames(df_r)
    data <- estimate_depth(df = data, j = 1:(nrow(data) - 1))
    data$h <- data$emax - data$emin
    return(data)
  }
  
  list_h <- lapply(X = x, FUN = estimate_height)
  
  df_h <- do.call(rbind, list_h)
  
  rownames(df_h) <- NULL
  
  return(df_h)
  
}

OCModel<-function (df, MA = 0, nwpath) {
  
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
    
    if (nrow(Pr)<4 | max(Pr$FAge)<MA) {
      next
    }
    
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

# global tendencies


#### Homogenize Age


TAll$FAge<-"NA"
TAll$FAge<-as.numeric(TAll$FAge)

for (i in 1:nrow(TAll)) {
  
  if (is.na(TAll[i,"Age"]) == FALSE) {TAll[i,"FAge"]<-TAll[i,"Age"]} 
  
  else {
    
    if (is.na(TAll[i,"Age.Pb"]) == FALSE) {TAll[i,"FAge"]<-TAll[i,"Age.Pb"]}
    
    else {if (is.na(TAll[i,"Age.C"]) == FALSE) {TAll[i,"FAge"]<-TAll[i,"Age.C"]}}}}


tendency_age<-as.data.frame(tendency(TAll, pnames="global")[1])
colnames(tendency_age)<-c("ID", "Ecosystem", "C_rho", "C_p", "C_Gr")

#### Count cores per group and get the percentages

NGr <- tendency_age %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

tendency_age %>% group_by(Ecosystem, C_Gr) %>% count()
NGrSg <- subset(tendency_age, Ecosystem == "Seagrass") %>% group_by(C_Gr) %>% count()
NGrSg %>% mutate(proc = ((n * 100) / sum(NGrSg[, 2])))

NGrSm <-
  subset(tendency_age, Ecosystem == "Tidal Marsh") %>% group_by(C_Gr) %>% count()
NGrSm %>% mutate(proc = ((n * 100) / sum(NGrSm[, 2])))

NGrMg <- subset(tendency_age, Ecosystem == "Mangrove") %>% group_by(C_Gr) %>% count()
NGrMg %>% mutate(proc = ((n * 100) / sum(NGrMg[, 2])))

# # TPb until 500 ---------------------------------------------------------


Folder = "Decay2023_Pb"
dir.create(Folder)

#### Homogenize Age


TPb$FAge<-"NA"
TPb$FAge<-as.numeric(TPb$FAge)

for (i in 1:nrow(TPb)) {
  
  if (is.na(TPb[i,"Age"]) == FALSE) {TPb[i,"FAge"]<-TPb[i,"Age"]} 
  
  else {TPb[i,"FAge"]<-TPb[i,"Age.Pb"]}}


  # model first 100 years ---------------------------------------------------
  
  
  #cores older than 80 years cut at 100
  
  Data_i<-subset(TPb, TPb$FAge < 100)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>80) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  # estimation of tendencies
  
  Data_t<-tendency(Data_i2, pnames="100")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_100Pb<-OCModel(DataAM, nwpath="Decay2023_Pb/100")
  
  
  
  
  
  
  
  # model 100-150 years ---------------------------------------------------
  
  Data_i<-subset(TPb, TPb$FAge < 150)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>100) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  
  Data_t<-tendency(Data_i2, pnames="100_150")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_150Pb<-OCModel(DataAM, MA= 100, nwpath="Decay2023_Pb/150")
  
  
  
  
  # model 150-300 years ---------------------------------------------------
  
  
  Data_i<-subset(TPb, TPb$FAge < 300)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>150) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  
  Data_t<-tendency(Data_i2, pnames="150_300")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_300Pb<-OCModel(DataAM, MA= 150, nwpath="Decay2023_Pb/300")
  
  
  # model 300-500 years ---------------------------------------------------
  
  Data_i<-subset(TPb, TPb$FAge < 500)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>300) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames="300_500")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_500Pb<-OCModel(DataAM, MA= 300, nwpath="Decay2023_Pb/500")
  
  
# # TC from 1000 to more than 2000 ---------------------------------------------------------


Folder = "Decay2023_C"
dir.create(Folder)

#### Homogenize Age


TC$FAge<-"NA"
TC$FAge<-as.numeric(TC$FAge)

for (i in 1:nrow(TC)) {
  
  if (is.na(TC[i,"Age"]) == FALSE) {TC[i,"FAge"]<-TC[i,"Age"]} 
  
  else {
    
    if (is.na(TC[i,"Age.C"]) == FALSE) {TC[i,"FAge"]<-TC[i,"Age.C"]}}}



  # model 500-1000 years ---------------------------------------------------
  
  
  Data_i<-subset(TC, TC$FAge < 1000)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>500) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames="500_1000")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_1000C<-OCModel(DataAM, MA= 500, nwpath="Decay2023_C/1000")
  
  
  
  
  # model 1000-1500 years ---------------------------------------------------
  
  
  Data_i<-subset(TC, TC$FAge < 1500)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>1000) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  
  Data_t<-tendency(Data_i2, pnames="1000_1500")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_1500C<-OCModel(DataAM, MA= 1000, nwpath="Decay2023_C/1500")
  
  
  
  # model 1500-2000 years ---------------------------------------------------
  
  
  Data_i<-subset(TC, TC$FAge < 2000)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>1500) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames="1500_2000")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_2000C<-OCModel(DataAM, MA= 1500, nwpath="Decay2023_C/2000")
  
  
  # model > 2000 years ---------------------------------------------------
  
  Data_i<-TC
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>2000) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames=">2000_C")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_m2000C<-OCModel(DataAM, MA= 2000, nwpath="Decay2023_C/more_2000")
  
  
# # TC and Pb from 1000 to more than 2000 ---------------------------------------------------------


Folder = "Decay2023_CPb"
dir.create(Folder)

#### Homogenize Age


TPbandC$FAge<-"NA"
TPbandC$FAge<-as.numeric(TPbandC$FAge)

for (i in 1:nrow(TPbandC)) {
  
  if (is.na(TPbandC[i,"Age"]) == FALSE) {TPbandC[i,"FAge"]<-TPbandC[i,"Age"]} }



  # model 500-1000 years ---------------------------------------------------
  
  
  Data_i<-subset(TPbandC, TPbandC$FAge < 1000)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>500) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames="500_1000")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_1000CPb<-OCModel(DataAM, MA= 500, nwpath="Decay2023_CPb/1000")
  
  
  
  
  # model 1000-1500 years ---------------------------------------------------
  
  
  Data_i<-subset(TPbandC, TPbandC$FAge < 1500)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>1000) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  
  Data_t<-tendency(Data_i2, pnames="1000_1500")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_1500CPb<-OCModel(DataAM, MA= 1000, nwpath="Decay2023_CPb/1500")
  
  
  
  # model 1500-2000 years ---------------------------------------------------
  
  
  Data_i<-subset(TPbandC, TPbandC$FAge < 2000)
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>1500) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames="1500_2000")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_2000CPb<-OCModel(DataAM, MA= 1500, nwpath="Decay2023_CPb/2000")
  
  
  # model > 2000 years ---------------------------------------------------
  
  Data_i<-TPbandC
  
  X<- split(Data_i, Data_i$Core)
  
  Data_i2<-data.frame(matrix(ncol = length(Data_i), nrow = 0))
  colnames(Data_i2)<-colnames(Data_i)
  
  for (i in 1:length(X)) {
    
    Data <- as.data.frame(X[i])
    colnames(Data)<-colnames(Data_i)
    
    if (max(Data$FAge)>2000) {
      
      Data_i2<-rbind(Data_i2, Data)
    }}
  
  
  Data_t<-tendency(Data_i2, pnames=">2000_C")
  
  TDEC<-Data_t[[2]]
  
  DataA <-
    as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])
  
  #carbon stock estimation por sample
  DataA<-estimate_h(DataA,
                    core = "Core",
                    mind = "Min.Depth",
                    maxd = "Max.Depth")
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
  
  #model
  DataAM<-as.data.frame(DataAM[, c("Core", "Ecosystem", "FAge", "Corg", "Corg.M")])
  
  fit_m2000CPb<-OCModel(DataAM, MA= 2000, nwpath="Decay2023_CPb/more_2000")
  
  
  
  
  
# eliminate models after visual check -------------------------------------

#eliminate some cores after visual check, we eliminate:

    #eliminate some cores after visual check
    fit_100Pb[c( 9, 10, 21, 23, 31), "k"]<-NA
    
    
    pairwise.wilcox.test(fit_100Pb$k, fit_100Pb$Ecosystem,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    fit_150Pb <- fit_150Pb[-c(6, 27, 33, 53), ]
    
    pairwise.wilcox.test(fit_150Pb$k, fit_150Pb$Ecosystem,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)    
    
    fit_300Pb <- fit_300Pb[-c(22, 23, 36), ]
    
    
    pairwise.wilcox.test(fit_300Pb$k, fit_300Pb$Ecosystem,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
    fit_500Pb <- fit_500Pb[-c(20, 21, 32), ]
    
    
    pairwise.wilcox.test(fit_500Pb$k, fit_500Pb$Ecosystem,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    
 
    #fit_1000C <- fit_1000C[-c( 5, 8, 16, 17, 18, 21, 28, 36, 46, 50, 51), ]
    
    #fit_1500C <- fit_1500C[-c(3, 10, 13, 19, 20, 21, 23, 30, 31, 32, 39, 43), ]
  
    #fit_2000C <- fit_2000C[-c( 2, 13, 14, 21, 22), ]

    #fit_m2000C <- fit_m2000C[-c(3, 5, 15, 16), ]
    #fit_m2000C <- fit_m2000C[!is.na(fit_m2000C$k),]
    
    # with C and Pb

    
    fit_1000CPb <- fit_1000CPb[-c( 5, 6, 10, 15, 21, 22, 23, 24), ]
    
    
    fit_1500CPb <- fit_1500CPb[-c(1, 2, 4, 5, 7, 9, 15, 16, 17), ]
    
    
    fit_2000CPb <- fit_2000CPb[-c( 1, 2, 6, 7), ]
    
    
    fit_m2000CPb <- fit_m2000C[-c(3), ]
    
    fit_m2000CPb <- fit_m2000CPb[!is.na(fit_m2000CPb$k),]  
    

# Final table and plots -------------------------------------------------------

Folder = "Decay2023"
    
    
k_table <-merge(fit_100Pb[,c(1,4)], fit_150Pb[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_300Pb[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_500Pb[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_1000CPb[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_1500CPb[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_2000CPb[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_m2000CPb[,c(1,4, 5)], by = "ID", all = TRUE)


colnames(k_table)<-c("ID", "k_100", "k_150", "k_300", "k_500", "k_1000","k_1500", "k_2000", "k_m2000", "Max_Age")
#k_table$k_100<-as.numeric(k_table$k_100)
#k_table$k_300<-as.numeric(k_table$k_300)

names(SingleCore)[names(SingleCore) == 'Core'] <- 'ID'

k_table<-merge(k_table, SingleCore[,c(1, 3, 7, 6, 8, 11, 12)], by = 'ID', all.x=T, all.y=F)
k_table<-merge(k_table, SAR[,c(1,7)], by = 'ID', all.x=T, all.y=F)

k_table_Mg<-subset(k_table, Ecosystem=='Mangrove')
k_table_Mg[,c(2:10)]<-sapply(k_table_Mg[,c(2:10)],FUN=as.numeric)
k_table_Sg<-subset(k_table, Ecosystem=='Seagrass')
k_table_Sg[,c(2:10)]<-sapply(k_table_Sg[,c(2:10)],FUN=as.numeric)
k_table_Sm<-subset(k_table, Ecosystem=='Tidal Marsh')
k_table_Sm[,c(2:10)]<-sapply(k_table_Sm[,c(2:10)],FUN=as.numeric)




### Table 1 ###


k_tablem<-melt (k_table[, c(1:8, 11)], id=c( "ID", "Ecosystem" ))




#normal distribution and significant differences among ecosystems

shapiro.test(k_table$k_150) #normal if pvalue > than 0.05



apply(k_table[,c(2:8)], FUN=shapiro.test, MARGIN = 2)


# differences among species

pairwise.wilcox.test(k_table$k_100, k_table$Ecosystem,
                     p.adjust.method = "BH")
temp<-subset(k_table, !Ecosystem=="Mangrove")
pairwise.wilcox.test(temp$k_1000, temp$Ecosystem,
                     p.adjust.method = "BH")

# summary table (manuscript Table 1) -----------------------------------



std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
cnt <- function(x) sum(!is.na(x))   

sum_table<-as.data.frame(colMeans(k_table_Mg[,c(2:8)], na.rm = TRUE))
sum_table[,2]<-as.data.frame(apply(k_table_Mg[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,3]<-as.data.frame(apply(k_table_Mg[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,4]<-colMeans(k_table_Sg[,c(2:8)], na.rm = TRUE)
sum_table[,5]<-as.data.frame(apply(k_table_Sg[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,6]<-as.data.frame(apply(k_table_Sg[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,7]<-colMeans(k_table_Sm[,c(2:8)], na.rm = TRUE)
sum_table[,8]<-as.data.frame(apply(k_table_Sm[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,9]<-as.data.frame(apply(k_table_Sm[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,10]<-colMeans(k_table[,c(2:8)], na.rm = TRUE)
sum_table[,11]<-as.data.frame(apply(k_table[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,12]<-as.data.frame(apply(k_table[,c(2:8)], FUN=cnt, MARGIN = 2))

colnames(sum_table)<-c("Mean Mangrove", "SE Mangrove", "n Mangrove","Mean Seagrass", "SE Seagrass", "n Seagrass",
                       "Mean Tidal Marsh", "SE Tidal Marsh", "n Tidal Marsh", "Mean All", "SE All", "n All")

Folder = "Decay2023"
write.csv(sum_table,
          file.path(Folder, "Summar decay.csv"),
          sep = ";",
          dec = ".")


# boxplot by timeframe figure ---------------------------------------------


k_tablem<-melt(k_table[,c(1:8, 11)], id = c("ID","Ecosystem")) 

k_tablem$variable <- as.character(k_tablem$variable)
k_tablem$variable[k_tablem$variable == 'k_100'] <- '80-100 yr'
k_tablem$variable[k_tablem$variable == 'k_150'] <- '100-150 yr'
k_tablem$variable[k_tablem$variable == 'k_300'] <- '150-300 yr'
k_tablem$variable[k_tablem$variable == 'k_500'] <- '300-500 yr'
k_tablem$variable[k_tablem$variable == 'k_1000'] <- '500-1000 yr'
k_tablem$variable[k_tablem$variable == 'k_1500'] <- '1000-1500 yr'
k_tablem$variable[k_tablem$variable == 'k_2000'] <- '1500-2000 yr'

k_tablem$value<-as.numeric(k_tablem$value)

ggplot(transform(k_tablem,
                 variable=factor(variable,levels=c('80-100 yr','100-150 yr','150-300 yr', '300-500 yr', '500-1000 yr', '1000-1500 yr', '1500-2000 yr'))),
       aes(Ecosystem, value))+ ggtitle("Decay rates by ecosystem and time frame")+ ylab("Decay rate (yr-1)") +
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))+
  facet_wrap(~variable)+
  scale_color_manual(values=c('blue', 'green4', "orange"))+
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0), 
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),)


k_tablem2<-subset(k_tablem, k_tablem$variable == "80-100 yr" | k_tablem$variable == "500-1000 yr")

box_100_1000<-  ggplot(transform(k_tablem2,
                                 variable=factor(variable,levels=c('80-100 yr',"500-1000 yr"))),
                       aes(Ecosystem, value))+ ggtitle("Decay rates by ecosystem and time frame")+ ylab("Decay rate (yr-1)") +
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))+
  facet_wrap(~variable)+
  scale_color_manual(values=c('blue', 'green4', "orange"))+
  ylim(0,0.0475)+
  theme(#legend.position = c(1, 0),
    #legend.justification = c(1, 0), 
    axis.title.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),)



ggsave(
  plot = box_100_1000,
  path = Folder,
  filename = "box_150_1000.jpg",
  units = "cm",
  width = 12,
  height = 7
)



# load decay rates from review --------------------------------------------

File <- "Data/k_rev.csv"

k_rev <- read.csv(File,
                  header = T,
                  sep = ";",
                  dec = ".")
k_rev <- as.data.frame(k_rev)
k_revS<- k_rev [-c(1,6,8,9),]
k_revM<- k_rev [c(8,9),]



# exponential model to predict k  -----------------------------------------

#fitting table

f_table_Sg<-as.data.frame(colMeans(subset(k_table, Ecosystem=="Seagrass")[,c(2:8)], na.rm = TRUE))
f_table_Sg[,2]<-c(90, 125, 225, 400, 750, 1250, 1750)
f_table_Sg[c(8:13),1]<-na.omit(k_table[,9])
f_table_Sg[c(8:13),2]<-na.omit(k_table[,10])
f_table_Sg[c(14:18), c(1,2)]<-k_revS[,c(2:3)]

#f_table_Sg<-f_table[-9,]
colnames(f_table_Sg)<-c("k","timeframe")

f_table_Tm<-as.data.frame(colMeans(subset(k_table, Ecosystem=="Tidal Marsh")[,c(2:6)], na.rm = TRUE))
f_table_Tm[,2]<-c(90, 125, 225, 400, 750)
f_table_Tm[c(6:7), c(1,2)]<-k_revM[,c(2:3)]
colnames(f_table_Tm)<-c("k","timeframe")

f_table_Mg<-as.data.frame(colMeans(subset(k_table, Ecosystem=="Mangrove")[,c(2:5)], na.rm = TRUE))
f_table_Mg[,2]<-c(50, 125, 225, 400)
colnames(f_table_Mg)<-c("k","timeframe")

# fit function k-timeframe

### exponential model to predict k 

kchange <- function(Tframe, A, C)
  (A * exp(C * Tframe))

      # seagrass model
      modelSg <-
        nls(
          k ~ kchange(timeframe, myA, myC),
          data = f_table_Sg,
          start = list(myA = 0.02, myC = 0.0003)
        )
      
      summary(modelSg)
      
      fitSg <- as.data.frame(c(1:5000))
      fitSg['new_col'] <- NA
      fitSg[, 2] <- kchange(c(1:5000), 0.04, -0.003)
      colnames(fitSg) <- list("timeframe", "predict")
      
      # tidal marsh model
      modelTm <-
        nls(
          k ~ kchange(timeframe, myA, myC),
          data = f_table_Tm,
          start = list(myA = 0.02, myC = 0.0003)
        )
      
      summary(modelTm)
      
      fitTm <- as.data.frame(c(1:5000))
      fitTm['new_col'] <- NA
      fitTm[, 2] <- kchange(c(1:5000), 0.02, -0.003)
      colnames(fitTm) <- list("timeframe", "predict")
      
      # Mangrove model
      modelMg <-
        nls(
          k ~ kchange(timeframe, myA, myC),
          data = f_table_Mg,
          start = list(myA = 0.02, myC = 0.0003)
        )
      
      summary(modelMg)
      
      fitMg <- as.data.frame(c(1:5000))
      fitMg['new_col'] <- NA
      fitMg[, 2] <- kchange(c(1:5000), 0.023, -0.0025)
      colnames(fitMg) <- list("timeframe", "predict")



# k vs time frame fitting figure (Figure 2) ------------------------------------------


std <- function(x) sd(x)/sqrt(length(x))


k_table_Mg<-subset(k_table, Ecosystem=='Mangrove')
k_table_Mg[,c(2:10)]<-sapply(k_table_Mg[,c(2:10)],FUN=as.numeric)
k_table_Sg<-subset(k_table, Ecosystem=='Seagrass')
k_table_Sg[,c(2:10)]<-sapply(k_table_Sg[,c(2:10)],FUN=as.numeric)
k_table_Sm<-subset(k_table, Ecosystem=='Tidal Marsh')
k_table_Sm[,c(2:10)]<-sapply(k_table_Sm[,c(2:10)],FUN=as.numeric)



fit_fig<-  
  ggplot(k_table_Sg, aes( Max_Age, k_m2000))+ ggtitle("Decay rate by time frame") + xlab("Time frame (years)") + ylab("Decay rate (yr-1)") +
  geom_point(color='green4')+
  geom_point(aes(50, mean(na.omit(k_table_Sg$k_100))), color='green4')+
  geom_errorbar(aes(50, ymin=mean(na.omit(k_table_Sg$k_100))-std(na.omit(k_table_Sg$k_100)), ymax=mean(na.omit(k_table_Sg$k_100))+std(na.omit(k_table_Sg$k_100))), color='green4')+
  geom_point(aes(50, mean(na.omit(k_table_Mg$k_100))), color='blue')+
  geom_errorbar(aes(50, ymin=mean(na.omit(k_table_Mg$k_100))-std(na.omit(k_table_Mg$k_100)), ymax=mean(na.omit(k_table_Mg$k_100))+std(na.omit(k_table_Mg$k_100))), color='blue')+
  geom_point(aes(50, mean(na.omit(k_table_Sm$k_100))), color='orange')+
  geom_errorbar(aes(50, ymin=mean(na.omit(k_table_Sm$k_100))-std(na.omit(k_table_Sm$k_100)), ymax=mean(na.omit(k_table_Sm$k_100))+std(na.omit(k_table_Sm$k_100))), color='orange')+
  
  geom_point(aes(125, mean(na.omit(k_table_Sg$k_150))), color='green4')+
  geom_errorbar(aes(125, ymin=mean(na.omit(k_table_Sg$k_150))-std(na.omit(k_table_Sg$k_150)), ymax=mean(na.omit(k_table_Sg$k_150))+std(na.omit(k_table_Sg$k_150))), color='green4')+
  geom_point(aes(125, mean(na.omit(k_table_Mg$k_150))), color='blue')+
  geom_errorbar(aes(125, ymin=mean(na.omit(k_table_Mg$k_150))-std(na.omit(k_table_Mg$k_150)), ymax=mean(na.omit(k_table_Mg$k_150))+std(na.omit(k_table_Mg$k_150))), color='blue')+  
  geom_point(aes(125, mean(na.omit(k_table_Sm$k_150))), color='orange')+
  geom_errorbar(aes(125, ymin=mean(na.omit(k_table_Sm$k_150))-std(na.omit(k_table_Sm$k_150)), ymax=mean(na.omit(k_table_Sm$k_150))+std(na.omit(k_table_Sm$k_150))), color='orange')+  
  
  geom_point(aes(175, mean(na.omit(k_table_Sg$k_300))), color='green4')+
  geom_errorbar(aes(175, ymin=mean(na.omit(k_table_Sg$k_300))-std(na.omit(k_table_Sg$k_300)), ymax=mean(na.omit(k_table_Sg$k_300))+std(na.omit(k_table_Sg$k_300))), color='green4')+
  geom_point(aes(175, mean(na.omit(k_table_Mg$k_300))), color='blue')+
  geom_errorbar(aes(175, ymin=mean(na.omit(k_table_Mg$k_300))-std(na.omit(k_table_Mg$k_300)), ymax=mean(na.omit(k_table_Mg$k_300))+std(na.omit(k_table_Mg$k_300))), color='blue')+
  geom_point(aes(175, mean(na.omit(k_table_Sm$k_300))), color='orange')+
  geom_errorbar(aes(175, ymin=mean(na.omit(k_table_Sm$k_300))-std(na.omit(k_table_Sm$k_300)), ymax=mean(na.omit(k_table_Sm$k_300))+std(na.omit(k_table_Sm$k_300))), color='orange')+
  
  geom_point(aes(400, mean(na.omit(k_table_Sg$k_500))), color='green4')+
  geom_errorbar(aes(400, ymin=mean(na.omit(k_table_Sg$k_500))-std(na.omit(k_table_Sg$k_500)), ymax=mean(na.omit(k_table_Sg$k_500))+std(na.omit(k_table_Sg$k_500))), color='green4')+
  geom_point(aes(400, mean(na.omit(k_table_Mg$k_500))), color='blue')+
  geom_errorbar(aes(400, ymin=mean(na.omit(k_table_Mg$k_500))-std(na.omit(k_table_Mg$k_500)), ymax=mean(na.omit(k_table_Mg$k_500))+std(na.omit(k_table_Mg$k_500))), color='blue')+
  geom_point(aes(400, mean(na.omit(k_table_Sm$k_500))), color='orange')+
  geom_errorbar(aes(400, ymin=mean(na.omit(k_table_Sm$k_500))-std(na.omit(k_table_Sm$k_500)), ymax=mean(na.omit(k_table_Sm$k_500))+std(na.omit(k_table_Sm$k_500))), color='orange')+
  
  geom_point(aes(750, mean(na.omit(k_table_Sg$k_1000))), color='green4')+
  geom_errorbar(aes(750, ymin=mean(na.omit(k_table_Sg$k_1000))-std(na.omit(k_table_Sg$k_1000)), ymax=mean(na.omit(k_table_Sg$k_1000))+std(na.omit(k_table_Sg$k_1000))), color='green4')+
  geom_point(aes(750, mean(na.omit(k_table_Mg$k_1000))), color='blue')+
  geom_errorbar(aes(750, ymin=mean(na.omit(k_table_Mg$k_1000))-std(na.omit(k_table_Mg$k_1000)), ymax=mean(na.omit(k_table_Mg$k_1000))+std(na.omit(k_table_Mg$k_1000))), color='blue')+
  geom_point(aes(750, mean(na.omit(k_table_Sm$k_1000))), color='orange')+
  geom_errorbar(aes(750, ymin=mean(na.omit(k_table_Sm$k_1000))-std(na.omit(k_table_Sm$k_1000)), ymax=mean(na.omit(k_table_Sm$k_1000))+std(na.omit(k_table_Sm$k_1000))), color='orange')+
  
  geom_point(aes(1250, mean(na.omit(k_table_Sg$k_1500))), color='green4')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(k_table_Sg$k_1500))-std(na.omit(k_table_Sg$k_1500)), ymax=mean(na.omit(k_table_Sg$k_1500))+std(na.omit(k_table_Sg$k_1500))), color='green4')+
  geom_point(aes(1250, mean(na.omit(k_table_Mg$k_1500))), color='blue')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(k_table_Mg$k_1500))-std(na.omit(k_table_Mg$k_1500)), ymax=mean(na.omit(k_table_Mg$k_1500))+std(na.omit(k_table_Mg$k_1500))), color='blue')+
  geom_point(aes(1250, mean(na.omit(k_table_Sm$k_1500))), color='orange')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(k_table_Sm$k_1500))-std(na.omit(k_table_Sm$k_1500)), ymax=mean(na.omit(k_table_Sm$k_1500))+std(na.omit(k_table_Sm$k_1500))), color='orange')+
  
  geom_point(aes(1750, mean(na.omit(k_table_Sg$k_2000))), color='green4')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(k_table_Sg$k_2000))-std(na.omit(k_table_Sg$k_2000)), ymax=mean(na.omit(k_table_Sg$k_2000))+std(na.omit(k_table_Sg$k_2000))), color='green4')+
  geom_point(aes(1750, mean(na.omit(k_table_Mg$k_2000))), color='blue')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(k_table_Mg$k_2000))-std(na.omit(k_table_Mg$k_2000)), ymax=mean(na.omit(k_table_Mg$k_2000))+std(na.omit(k_table_Mg$k_2000))), color='blue')+
  geom_point(aes(1750, mean(na.omit(k_table_Sm$k_2000))), color='orange')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(k_table_Sm$k_2000))-std(na.omit(k_table_Sm$k_2000)), ymax=mean(na.omit(k_table_Sm$k_2000))+std(na.omit(k_table_Sm$k_2000))), color='orange')+
  
  geom_point(data= k_revS,mapping = aes(k_revS$Max_Age, k_revS$k), color='green', shape= 17)+
  geom_point(data= k_revM,mapping = aes(k_revM$Max_Age, k_revM$k), color='orange4', shape= 17)+
  
  geom_line(data= fitSg, mapping = aes(timeframe, predict), color='green4')+
  geom_line(data= fitTm, mapping = aes(timeframe, predict), color='orange')+
  geom_line(data= fitMg, mapping = aes(timeframe, predict), color='blue')+
  
  xlim(0,4000)+
  
  annotate("text", x=1500, y=0.025, color= "blue",  size = 5, label= expression(y == 0.03 * e ** (-0.0025 * 
                                                                                                     x)))+
  annotate("text", x=1500, y=0.02, color= "green4",size = 5,label= expression(y == 0.04 * e ** (-0.003 * 
                                                                                                   x)))+
  annotate("text", x=1500, y=0.015, color= "orange",size = 5,label= expression(y == 0.02 * e ** (-0.003 * 
                                                                                                    x)))


ggsave(
  plot = fit_fig,
  path = Folder,
  filename = "fit_plot.jpg",
  units = "cm",
  width = 12,
  height = 7
)



# correlation decay rate, OC, Mud and d13C (Figure 3) ------------------------------------

  # first we estimate the average OC, mud, d13C and SAR content for the studied time frame
  # we use the max age of each core
  
  
  #write a function to estimate correlation between function and k for specific time frames
  
  
  
  estimate_sum_var <-function (df, df2) {
    
    x<-split(df, df$Core)
    
    temp<-df2
    temp$OC<-"NA"
    temp$Mud<-"NA"
    temp$d13C<-"NA"
    temp$SAR<-"NA"
    
    
    for (i in 1:nrow(temp)) {
      
      data<-as.data.frame(x[temp[i,"ID"]])
      colnames(data)<-colnames(df)
      data_a<-subset(data, data$FAge<=temp[i, "Max.Age"])
      temp[i, "OC"]<-mean(data_a$Corg, na.rm=TRUE)
      temp[i, "Mud"]<-mean(data_a$Mud, na.rm=TRUE)
      temp[i, "d13C"]<-mean(data_a$d13C, na.rm=TRUE)
      temp[i, "SAR"]<-max(data_a$Max.Depth)/max(data_a$FAge)}
    
    temp$Mud<-as.numeric(temp$Mud)
    temp$OC<-as.numeric(temp$OC)
    temp$d13C<-as.numeric(temp$d13C)
    temp$SAR<-as.numeric(temp$SAR)
    
    return(temp)
    
    
  }
  



  var_100<-estimate_sum_var (TPb, fit_100Pb)
  var_100 <-merge(var_100, oracle[,c(1,3,4,7:12)], by = "ID", all = T)
  var_100 <- var_100[!is.na(var_100$k), ]
  
  write.csv(var_100,file.path(Folder,"var_100.csv"),sep=";", dec=",")
  
  
  var_1000<-estimate_sum_var (TPbandC, fit_1000CPb)
  var_1000 <-merge(var_1000, oracle[,c(1,3,4,7:12)], by = "ID", all = T)
  var_1000 <- var_1000[!is.na(var_1000$k), ]
  
  write.csv(var_1000,file.path(Folder,"var_1000.csv"),sep=";", dec=",")
  

# PCA decay ---------------------------------------------------------------

  temp<-var_1000[,c(2,6,8,9,11:17)]
  
  temp<-na.omit(temp)
  
  
  
  pca<-as.data.frame(scale(temp[,-1]))
  
  
  # PCA with varimax rotation #
  PCA<-principal(pca, nfactors=1, residual=T, rotate="varimax", covar=F) #Modify number of factor until you have the maximun number of components with more than 1 explained variance
  AutVal<-sum(PCA$values>1)
  PCA<-principal(pca, nfactors=AutVal, residual=T, rotate="varimax", covar=F)
  print(PCA)
  
  fa.parallel(pca)
  
  
  # loading #
  loa <-loadings(PCA)
  loa <-as.data.frame(loa[,1:ncol(loa)])
  
  # comunality #
  
  loa2 <- loa^2 #Cuadrados de los loadings
  Com <-t(loa2) # Transpuesta
  
  barplot(Com, col=c("skyblue1","tan1","springgreen2","lightgoldenrod2","violetred3","aquamarine", "chartreuse4"),
          ylim=c(0,1), axisnames = T, las=2, cex.names=0.8,xaxt="s",main="Fraccionamiento comunalidad")# barplot
  
  write.csv(loa,file.path(Folder,"1000_short_loa.csv"),sep=";", dec=",")
  
  
  
  # scores #
  sco<-PCA$scores
  sco<-as.data.frame(sco)
  write.csv(sco,file.path(Folder,"1000_short_sco.csv"), sep=";", dec=".")
  
  
  ########### biplot carmen towapoh pa clr #################
  
  colnames(sco)<-c("PC1","PC2","PC3","PC4") ###"Depth"
  colnames(loa)<-c("PC1","PC2","PC3", "PC4")
  fit<- as.data.frame(sco)
  Eco<-temp[,1]
  
  
  col_vector <- c("blue","orange", "green4")
  
  ggplot(data=loa, aes(x=PC1,y=PC3)) +geom_point(data=fit,aes(x=PC1,y=PC3,colour=Eco), size=2)+ coord_fixed(ratio=1)+ theme_bw() +
    annotate("text", x=(loa$PC1*5.8), y=(loa$PC3*5.8), label=row.names(loa), col="red3", size= 5)+
    xlab("PC1 (22%)")+ylab("PC2 (16%)")+ #theme(legend.position = "none")+
    #xlim(c(-4.8,4.8))+
    # ylim(c(-3,5))+
    scale_color_manual(values = col_vector)+
    scale_shape_manual(values=c(15,17,19))+
    #scale_shape_manual(values=c(20))+
    theme(legend.title = element_text(colour = "white"),
          legend.position = "bottom")+
    geom_segment(data=loa, aes(x =0 , y = 0, xend = PC1*5, yend = PC3*5), arrow = arrow(length = unit(1/2, 'picas')), color = "black", linetype='solid', size=0.5)
  
  #ggsave("_Gr3_1-5.jpg", units="cm", width = 20, height = 20)
  
  
  # map of long_100 and short_1000 decays
  
  ggplot() + xlab("Longitude") + ylab("Latitude") +
    geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
    #geom_point(aes(x = long, y = lat))+
    geom_point(data= var_100,aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
    coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
    scale_fill_manual(values = c("blue",  "green","orange")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot() + xlab("Longitude") + ylab("Latitude") +
    geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
    #geom_point(aes(x = long, y = lat))+
    geom_point(data= var_1000,aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
    coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
    scale_fill_manual(values = c(  "green","orange")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  

# map of fitted cores (Figure 5) ----------------------------------------------------
  
  
  k_table_c<-k_table[rowSums(is.na(k_table[,c(2:9)])) != ncol(k_table[,c(2:9)]), ]
  
  
  fitted_cores<-SingleCore[SingleCore$ID %in% k_table_c$ID,]
  
  
  #load a world map
  WM <- map_data("world")
  
  global<- fitted_cores %>%
    ggplot() + ggtitle("Sampling sites for OC decay modeled cores") + xlab("Longitude") + ylab("Latitude") +
    geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
    #geom_point(aes(x = long, y = lat))+
    geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
    coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
    scale_fill_manual(values = c("blue",  "green","orange")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  nam<-fitted_cores %>%
    ggplot() + xlab("Longitude") + ylab("Latitude") +
    geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
    #geom_point(aes(x = long, y = lat))+
    geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
    coord_sf(xlim = c(-150, -50), ylim = c(-20, 80)) +
    scale_fill_manual(values = c("blue",  "green","orange")) +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none")
  
  
  eu<-fitted_cores %>%
    ggplot()  + xlab("Longitude") + ylab("Latitude") +
    geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
    #geom_point(aes(x = long, y = lat))+
    geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
    coord_sf(xlim = c(-10, 50), ylim = c(20, 60)) +
    scale_fill_manual(values = c("blue",  "green","orange")) +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none")
  
  
  aus<-fitted_cores %>%
    ggplot()  + xlab("Longitude") + ylab("Latitude") +
    geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
    #geom_point(aes(x = long, y = lat))+
    geom_point(aes(x = Long, y = Lat,  fill = Ecosystem), pch = 21, size = 1.8) +
    coord_sf(xlim = c(110, 155), ylim = c(-40, -5))  +
    scale_fill_manual(values = c("blue",  "green","orange")) +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none")
  
  ssfc<-grid.arrange(global, nam, eu, aus, 
                    layout_matrix = rbind(c(1, 1, 1),
                                          c(2, 3, 4)))
  
  
  ggsave(
    plot = ssfc,
    path = Folder,
    filename =  "Sampling sites fitted cores.jpg",
    units = "cm",
    width = 20,
    height = 15
  )
  
  

# decay rates first meter of the soil -------------------------------------------------------------
  

  
  extract_depth_by_age<- function (df, AGE = 100) {
    
    core <- df[1,"Core"]
    c_length <- max(df[,"Max.Depth"])
    age_depth <- df[which.min(abs(AGE-df$Max.Depth)),"FAge"]
    
    output <- data.frame(core = core, Length= c_length, age_depth = age_depth)
    
    return(output)
    
  }
  
  x<-split(TAll, TAll$Core)
  
  core_100_age_l <- lapply(X = x,  extract_depth_by_age, AGE = 100) # return a list
  core_100_age <- as.data.frame(do.call(rbind, core_100_age_l)) # from list to dataframe
  
  core_100_age <- subset(core_100_age, core_100_age$Length > 90)
  
  mean(core_100_age$age_depth)
  std(core_100_age$age_depth)
  
  
  
  # from time to decay rate
  
  # degradation of the first meter
  #seagrass
  0.04 * exp(-0.003 * 1713.38)
  
  0.04 * exp(-0.003 * 1401.71)
  0.04 * exp(-0.003 * 2025.05)
  
  #mangrove
  0.028 * exp(-0.004 * 1713.38)
  
  0.028 * exp(-0.004 * 1401.71)
  0.028 * exp(-0.004 * 2025.05)
  
  #tidal marshes
  0.026 * exp(-0.0025 * 1713.38)
  
  0.026 * exp(-0.0025 * 1401.71)
  0.026 * exp(-0.0025 * 2025.05)
  
  
  # depth of the first 100 and 1000 years -------------------------------------------------------------
  
  
  
  extract_age_by_depth<- function (df, DEPTH = 100) {
    
    core <- df[1,"Core"]
    c_age <- max(df[,"FAge"])
    depth_age <- df[which.min(abs(DEPTH-df$FAge)),"Max.Depth"]
    
    output <- data.frame(core = core, Age= c_age, depth_age = depth_age)
    
    return(output)
    
  }
  
  x<-split(TAll, TAll$Core)
  
  core_100_depth_l <- lapply(X = x,  extract_age_by_depth, DEPTH = 100) # return a list
  core_100_depth <- as.data.frame(do.call(rbind, core_100_depth_l)) # from list to dataframe
  
  mean(core_100_depth$depth_age)
  std(core_100_depth$depth_age)
  
  
  
  
  

  

