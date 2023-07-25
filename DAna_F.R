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

############## Loading the dataset #####################

File <- "Data/Cores.csv"

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

B <- subset(Cores, Cores$V.vs.B == "Vegetated" | Cores$V.vs.B == "vegetated")
unique(B[, 5])
#Bare <- subset(A, A[,5]=='Bare')
length(unique(B$Core))
SingleCore<-B[!duplicated(B$Core),]

table(SingleCore$Ecosystem)


B$Corg <- as.numeric(B$Corg)
B$Mud <- as.numeric(B$Mud)
B$DBD <- as.numeric(B$DBD)

###### Sampling sites Map #########

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


ggplot(CM, aes(Ecosystem, Av_Mud_25)) +
  geom_boxplot() +
  geom_jitter()

ggplot(CM, aes(Ecosystem, Av_13C_25)) +
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
CM_nMg<-subset(CM, !Ecosystem=="Mangrove")
pairwise.wilcox.test(CM_nMg$Av_13C_25, CM_nMg$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

write.csv(CM,
          file.path(Folder, "AvMdC.csv"),
          sep = ";",
          dec = ".")



# Sediment accretion rate -------------------------------------------------

File <- "Data/Acc_Mass-Age_F.csv"

Dates <- read.csv(File,
                  header = T,
                  sep = ",",
                  dec = ".")
Dates <- as.data.frame(Dates)

# estimates SAR last 150 years 

SAR <- data.frame(
  ID = character(),
  Ecosystem = character(),
  SAR_Age = numeric(),
  SAR_Pb = numeric(),
  SAR_C = numeric())

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
  
}

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

##### Av Corg and trends

ggplot(CM, aes(DT$C_Gr, CM$Av_C)) +
  geom_boxplot() +
  geom_jitter()

names(SingleCore)[names(SingleCore) == 'Core'] <- 'ID'
DT<-merge(DT, SingleCore[,c(1, 3, 7, 6, 8, 11, 12)], by = 'ID', all.x=T, all.y=F)

library(janitor)
counts<-tabyl(DT, Specie, C_Gr)
counts<-counts %>% mutate(percent = counts[,2]*100/((counts[,2])+(counts[,3])+(counts[,4])))

counts[c(2,8,13,14,15,18:22),6]<-"Persistent"
counts[c(4:7, 17, 23:25),6]<-"Opportunistic"
counts[c(9:12, 16),6]<-"Colonising"

ggplot(counts,aes(V6, percent))+
  geom_boxplot()+
  geom_jitter()


counts<-tabyl(DT, Life.form, C_Gr)
counts<-counts %>% mutate(percent = counts[,2]*100/((counts[,2])+(counts[,3])+(counts[,4])))

#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

shapiro.test(CM$Av_C_25) #(>0.05 normal, <0.05 no normal)

## Student's t-test  if normally distributed, wilcox if not


DT2<-DT
DT2$C_Gr <- recode(DT2$C_Gr, DEC = 'DEC',
                  INC  = 'N',
                  NT = 'N')

ggplot(CM, aes(DT2$C_Gr, CM$Av_C_25)) +
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(CM$Av_C_25, DT2$C_Gr,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

DT2<-cbind(CM, DT2)
DT2<-DT2[,c(1:20, 25)]

#Seagrass meadows
    DT2Sg<-subset(DT2, Ecosystem=="Seagrass")
    #organic carbon  
    
    SC<-ggplot(DT2Sg, aes(C_Gr, Av_C_25)) + ggtitle("Seagrass")+
        geom_boxplot()+
        geom_jitter(color="green4")+
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank())
              
      
      pairwise.wilcox.test(DT2Sg$Av_C_25, DT2Sg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)
    #mud  
      SM<-ggplot(DT2Sg, aes(C_Gr, Av_Mud_25)) + 
        geom_boxplot()+
        geom_jitter(color="green4")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
        
      
      pairwise.wilcox.test(DT2Sg$Av_Mud_25, DT2Sg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)  
      
     
      
      
       
#Tidal marshes      
      #organic carbon  
      TMC<-ggplot(DT2Sm, aes(C_Gr, Av_C_25)) + ggtitle("Tidal Marsh")+
        geom_boxplot()+
        geom_jitter(color="orange")+
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      pairwise.wilcox.test(DT2Sm$Av_C_25, DT2Sm$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)
      #mud  
      TMM<-ggplot(DT2Sm, aes(C_Gr, Av_Mud_25)) +
        geom_boxplot()+
        geom_jitter(color="orange")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      pairwise.wilcox.test(DT2Sm$Av_Mud_25, DT2Sm$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)  
      
#Mangroves        
    DT2Mg<-subset(DT2, Ecosystem=="Mangrove")
    #organic carbon  
    MGC<-ggplot(DT2Mg, aes(C_Gr, Av_C_25)) + ylab("OC% (Top 25cm)") + ggtitle("Mangrove")+
      geom_boxplot()+
      geom_jitter(color="blue")+
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank())
    
    pairwise.wilcox.test(DT2Mg$Av_C_25, DT2Mg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05)
    #mud  
    TGM<-ggplot(DT2Mg, aes(C_Gr, Av_Mud_25)) + ylab("Mud% (Top 25cm)") +
      geom_boxplot()+
      geom_jitter(color="blue")+
      theme(axis.title.x = element_blank())
    
    pairwise.wilcox.test(DT2Mg$Av_Mud_25, DT2Mg$C_Gr,
                         p.adjust.method = "BH") # are significantly different (p < 0.05) 


    
todos_CM<-grid.arrange(MGC, SC, TMC, TGM, SM, TMM, ncol=3)
    
ggsave( plot = todos_CM,
      path = Folder,
       filename =  "C_Mud_Gr.jpg",
       units = "cm",
       width = 20,
       height = 10
)    




#### Count cores per group and get the percentages

NGr <- DT %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

DT %>% group_by(Ecosystem.x, C_Gr) %>% count()
NGrSg <- subset(DT, Ecosystem.x == "Seagrass") %>% group_by(C_Gr) %>% count()
NGrSg %>% mutate(proc = ((n * 100) / sum(NGrSg[, 2])))

NGrSm <-
  subset(DT, Ecosystem.x == "Tidal Marsh") %>% group_by(C_Gr) %>% count()
NGrSm %>% mutate(proc = ((n * 100) / sum(NGrSm[, 2])))

NGrMg <- subset(DT, Ecosystem.x == "Mangrove") %>% group_by(C_Gr) %>% count()
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

File <- "Data/Acc_Mass-Age_F.csv"

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



TAll = filter(C, !is.na(Age) | !is.na(Age.Pb) | !is.na(Age.C))# df with core with any model
length(unique(TAll$Core))

TPb = subset(C, !is.na(Age) | !is.na(Age.Pb) )# df with core with Pb model
length(unique(TPb$Core))

TC = filter(C, !is.na(Age.C))# df with core with C model
length(unique(TC$Core))


TPbandC = filter(C, !is.na(Age) | !is.na(Age.Pb))# df with core with C and Pb model
TPbandC = filter(TPbandC, !is.na(Age.C))
length(unique(TPbandC$Core))


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
    filename = paste(pnames,'_TNT.jpg'),
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
    filename = paste(pnames,'_TDEC.jpg'),
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
    filename = paste(pnames,'_TINC.jpg'),
    width = 20,
    height = 30,
    units = 'cm'
  )
  
  df_list<-list(TDT, TDEC, TINC, TNT)
  names(df_list)<-c("Sp_TDT", "TDEC", "TINC", "TNT")
  return(df_list)
}
  
prueba2<-tendency(TAll, pnames="prueba")

TDEC<-prueba2[[2]]

tendAll<-prueba2[[1]]

#### Count cores per group and get the percentages

NGr <- tendAll %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

tendAll %>% group_by(Ecosystem, C_Gr) %>% count()
NGrSg <- subset(tendAll, Ecosystem == "Seagrass") %>% group_by(C_Gr) %>% count()
NGrSg %>% mutate(proc = ((n * 100) / sum(NGrSg[, 2])))

NGrSm <-
  subset(tendAll, Ecosystem == "Tidal Marsh") %>% group_by(C_Gr) %>% count()
NGrSm %>% mutate(proc = ((n * 100) / sum(NGrSm[, 2])))

NGrMg <- subset(tendAll, Ecosystem == "Mangrove") %>% group_by(C_Gr) %>% count()
NGrMg %>% mutate(proc = ((n * 100) / sum(NGrMg[, 2])))




### groups and SAR ###

#Prueba<-left_join(TDT,SAR, "ID")

#ggplot(Prueba, aes(C_Gr, SAR))+
#  geom_boxplot()+
 # geom_jitter(aes(color=Ecosystem))

#ggplot(Prueba, aes(C_Gr, SAR))+
 # geom_boxplot(aes(color=Ecosystem))+
  #geom_jitter(aes(color=Ecosystem))


#pairwise.wilcox.test(Prueba$SAR, Prueba$C_Gr,
 #                    p.adjust.method = "BH") # are significantly different (p < 0.05)


#################################
### NLM time-Accumulated Mass ###
#################################


#For those cores that decrease with time

DataA <-
  as.data.frame(TDEC[, c("Core", "Ecosystem", "DBD","Min.Depth","Max.Depth","FAge", "Corg")])


#Estimate organic carbon accumulated mass

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

DataA<-estimate_h(DataA,
                  core = "Core",
                  mind = "Min.Depth",
                  maxd = "Max.Depth")


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

# to fit the model the core has to have more than 3 points
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

max_depth<-OCModel(DataAM, nwpath="Decay2023/Max_Depth")

#eliminate outlayers:
#eliminate empty cores (no model): 
#eliminate some cores after visual check, we eliminate: 
max_depth <- max_depth[-c(), ]



ggplot(TfitsM_DEC, aes(x = k)) +
  geom_histogram()

shapiro.test(TfitsM_DEC$k) #(>0.05 normal, <0.05 no normal)

# check if time frame has an effect 

ggplot(TfitsM_DEC, aes(Max.Age, k))+
  geom_point()



# model first 100 years ---------------------------------------------------

Data_i<-subset(TAll, TAll$FAge < 100)
Data_t<-tendency(Data_i, pnames="100")

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

fit_100<-OCModel(DataAM, nwpath="Decay2023/100")

#eliminate some cores after visual check, we eliminate: Sg_081, Sg_111, Sg_241, Sg_316, Sg_317, Sg_321, Sg_323, Sg_332, 
#Sg_497, Sm_004, Sm_68, Sm_069, Sm_092, Sm_97, Sm105
fit_100[c(13, 17, 27, 38, 39, 40, 41, 42, 47, 49, 61, 62, 67, 69, 70), "k"]<-NA


ggplot(fit_100, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

ggplot(fit_100, aes( Max.Age, k))+
  geom_point()


pairwise.wilcox.test(fit_100$k, fit_100$Ecosystem,
                    p.adjust.method = "BH") # are significantly different (p < 0.05)



# model 100-150 years ---------------------------------------------------

Data_i<-subset(TAll, TAll$FAge < 150)
Data_t<-tendency(Data_i, pnames="100_150")

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

fit_150<-OCModel(DataAM, MA= 100, nwpath="Decay2023/150")


#eliminate outlayers: 

ggplot(fit_150, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Mg_023, Sg_193, Sm_004, Sm_049
fit_150 <- fit_150[-c(6, 32, 59, 63), ]

#eliminate empty cores (no model)
fit_150<-fit_150[!is.na(fit_100$P),]


ggplot(fit_150, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(fit_150$k, fit_150$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model 150-300 years ---------------------------------------------------


Data_i<-subset(TAll, TAll$FAge < 300)
Data_t<-tendency(Data_i, pnames="150_300")

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

fit_300<-OCModel(DataAM, MA= 150, nwpath="Decay2023/300")


#eliminate outlayers: 

ggplot(fit_300, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_179, Sg_192, Sg_194, Sg_250, Sg_476, Sg_485, Sg_495, Sg_497
#Sm_004
fit_300 <- fit_300[-c( 38, 40, 41, 42, 46, 61, 68, 70, 72, 73), ]

#eliminate empty cores (no model)
fit_300<-fit_300[!is.na(fit_300$P),]


ggplot(fit_300, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()


pairwise.wilcox.test(fit_300$k, fit_300$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


# model 300-500 years ---------------------------------------------------

Data_i<-subset(TAll, TAll$FAge < 500)
Data_t<-tendency(Data_i, pnames="300_500")

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

fit_500<-OCModel(DataAM, MA= 300, nwpath="Decay2023/500")


#eliminate outlayers: 

ggplot(fit_500, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_041, Sg_193, Sg_250, Sg_310, Sg_312, Sg_314, Sg_476, Sg_479, 
#Sg_480, Sg_495, Sg_497, Sm_004
fit_500 <- fit_500[-c(15, 42, 47, 51, 53, 54, 61, 64, 65, 73, 75, 76 ), ]

#eliminate empty cores (no model)
fit_500<-fit_500[!is.na(fit_500$P),]


ggplot(fit_500, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(fit_500$k, fit_500$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model 500-1000 years ---------------------------------------------------


Data_i<-subset(TAll, TAll$FAge < 1000)
Data_t<-tendency(Data_i, pnames="500_1000")

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

fit_1000<-OCModel(DataAM, MA= 500, nwpath="Decay2023/1000")


#eliminate outlayers: 

ggplot(fit_1000, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_019, Sg_041, Sg_111, Sg_117, Sg_195, Sg_279, Sg_312, Sg_314, 
#Sg_490, Sm_004, Sm_010
fit_1000 <- fit_1000[-c( 12, 15, 27, 32, 45, 50, 55, 56, 73, 80, 81), ]

#eliminate empty cores (no model)
fit_1000<-fit_1000[!is.na(fit_1000$P),]


ggplot(fit_1000, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()


pairwise.wilcox.test(fit_1000$k, fit_1000$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model 1000-1500 years ---------------------------------------------------


Data_i<-subset(TAll, TAll$FAge < 1500)
Data_t<-tendency(Data_i, pnames="1000_1500")

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

fit_1500<-OCModel(DataAM, MA= 1000, nwpath="Decay2023/1500")


#eliminate outlayers: 

ggplot(fit_1500, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_041, Sg_097, Sg_112, Sg_121, Sg_170, Sg_195, Sg_311, Sg_312, 
#Sg_314, Sg_490, Sm_004, Sm_022
fit_1500 <- fit_1500[-c( 15, 25, 28, 35, 37, 45, 54, 55, 56, 73, 81, 83), ]

#eliminate empty cores (no model)
fit_1500<-fit_1500[!is.na(fit_1500$P),]


ggplot(fit_1500, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(fit_1500$k, fit_1500$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


# model 1500-2000 years ---------------------------------------------------


Data_i<-subset(TAll, TAll$FAge < 2000)
Data_t<-tendency(Data_i, pnames="1500_2000")

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

fit_2000<-OCModel(DataAM, MA= 1500, nwpath="Decay2023/2000")


#eliminate outlayers: 

ggplot(fit_2000, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_041, Sg_170, Sg_191, Sg_312, Sg_314
fit_2000 <- fit_2000[-c( 15, 35, 40, 54, 55), ]

#eliminate empty cores (no model)
fit_2000<-fit_2000[!is.na(fit_2000$P),]


ggplot(fit_2000, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()


pairwise.wilcox.test(fit_2000$k, fit_2000$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model > 2000 years ---------------------------------------------------


Data_i<-TAll
Data_t<-tendency(Data_i, pnames=">2000")

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

fit_m2000<-OCModel(DataAM, MA= 2000, nwpath="Decay2023/more_2000")


#eliminate outlayers: 

ggplot(fit_m2000, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_041, Sg_097, Sg_314
fit_m2000 <- fit_m2000[-c( 17, 27, 54), ]

#eliminate empty cores (no model)
fit_m2000<-fit_m2000[!is.na(fit_m2000$P),]


ggplot(fit_m2000, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()



# Final table and plots -------------------------------------------------------

k_table <-merge(fit_100[,c(1,4)], fit_150[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_300[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_500[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_1000[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_1500[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_2000[,c(1,4)], by = "ID", all = TRUE)
k_table <-merge(k_table, fit_m2000[,c(1,4, 5)], by = "ID", all = TRUE)

colnames(k_table)<-c("ID", "k_100", "k_150", "k_300", "k_500", "k_1000","k_1500", "k_2000", "k_m2000", "Max_Age")
k_table$k_100<-as.numeric(k_table$k_100)
k_table$k_300<-as.numeric(k_table$k_300)

names(SingleCore)[names(SingleCore) == 'Core'] <- 'ID'

k_table<-merge(k_table, SingleCore[,c(1, 3, 7, 6, 8, 11, 12)], by = 'ID', all.x=T, all.y=F)
k_table<-merge(k_table, SAR[,c(1,5)], by = 'ID', all.x=T, all.y=F)


ggplot(temp, aes(k_150, Bioregions))+
  geom_boxplot()+
  geom_jitter(aes(color=temp$Specie))


ggplot(temp, aes(k_100, Specie))+
  geom_boxplot()+
  geom_jitter(aes(color=temp$Bioregions))


#normal distribution and significant differences among ecosystems

shapiro.test(k_table$k_100) #normal if pvalue > than 0.05

apply(k_table[,c(2:8)], FUN=shapiro.test, MARGIN = 2)


temp<-subset(k_table_Sg, !is.na(k_100))

pairwise.wilcox.test(temp$k_100, temp$Bioregions,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


# summary table (manuscript Table 1)

k_table_Mg<-subset(k_table, Ecosystem=='Mangrove')
k_table_Mg<-sapply(k_table_Mg,FUN=as.numeric)
k_table_Sg<-subset(k_table, Ecosystem=='Seagrass')
k_table_Sg<-sapply(k_table_Sg,FUN=as.numeric)
k_table_Sm<-subset(k_table, Ecosystem=='Tidal Marsh')
k_table_Sm<-sapply(k_table_Sm,FUN=as.numeric)

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


write.csv(sum_table,
          file.path(Folder, "Summar decay.csv"),
          sep = ";",
          dec = ".")



# load decay rates from review

File <- "Data/k_rev.csv"

k_rev <- read.csv(File,
                  header = T,
                  sep = ";",
                  dec = ".")
k_rev <- as.data.frame(k_rev)
k_rev<- k_rev [-c(1,6),]

# boxplot by timeframe figure

mk_table<-melt(k_table[,-c(9,10)], id = c("ID","Ecosystem")) 

mk_table$variable <- as.character(mk_table$variable)
mk_table$variable[mk_table$variable == 'k_100'] <- '0-100 yr'
mk_table$variable[mk_table$variable == 'k_150'] <- '100-150 yr'
mk_table$variable[mk_table$variable == 'k_300'] <- '150-300 yr'
mk_table$variable[mk_table$variable == 'k_500'] <- '300-500 yr'
mk_table$variable[mk_table$variable == 'k_1000'] <- '500-1000 yr'
mk_table$variable[mk_table$variable == 'k_1500'] <- '1000-1500 yr'
mk_table$variable[mk_table$variable == 'k_2000'] <- '1500-2000 yr'

ggplot(transform(mk_table,
                 variable=factor(variable,levels=c('0-100 yr','100-150 yr','150-300 yr', '300-500 yr', '500-1000 yr', '1000-1500 yr', '1500-2000 yr'))),
       aes( Ecosystem, value))+ ggtitle("Decay rates by ecosystem and time frame")+ ylab("Decay rate (yr-1)") +
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))+
  facet_wrap(~variable)+
  scale_color_manual(values=c('blue', 'green4', "orange"))+
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0), 
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),)


#fitting table

f_table<-as.data.frame(colMeans(k_table[,c(2:8)], na.rm = TRUE))
f_table[,2]<-c(100, 150, 300, 500, 1000, 1500, 2000)
f_table[c(8:12),1]<-na.omit(k_table[,9])
f_table[c(8:12),2]<-na.omit(k_table[,10])
f_table[c(13:16), c(1,2)]<-k_rev[,c(2:3)]


# fit function k-timeframe

### exponential model to predict k in Posidonia meadows


kchange <- function(Tframe, A, C)
  (A * exp(C * Tframe))

model1 <-
  nls(
    f_table[,1] ~ kchange(f_table[,2], myA, myC),
    data = f_table,
    start = list(myA = 0.04, myC = -0.00201)
  )


plot(f_table[,1], f_table[,2])

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




# correlation decay rate (150) and SAR ------------------------------------


plot(k_table_Sg$k_100, k_table_Sg$SAR)

ggplot(k_table_Sg,aes(k_100, SAR))+
  geom_point(aes(color=Bioregions))

cor.test(k_table$k_100, k_table$SAR, method=c("pearson"))

#check for outlayers (https://www.r-bloggers.com/2016/12/outlier-detection-and-treatment-with-r/)

# visual check
boxplot(k_100 ~ Ecosystem, data=k_table, main="Decay 150yr by Ecosystem")  # clear pattern is noticeable.
boxplot(SAR ~ Ecosystem, data=k_table, main="SAR by Ecosystem")  # this may not be significant, as day of week variable is a subset of the month var.

#Cook’s Distance
mod <- lm(SAR ~ k_100, data=k_table)
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

# one outlayer: 100

cor.test(k_table[,-c(51,100)]$k_100, k_table[,-c(51,100)]$SAR, method=c("pearson"))

ggplot(k_table[,-c(51,100)],aes(k_100, SAR))+
  geom_point(aes(color=Ecosystem))





















# k vs time frame fitting figure

std <- function(x) sd(x)/sqrt(length(x))

ggplot(k_table_Sg, aes( Max_Age, k_m2000))+ ggtitle("Decay rate by time frame") + xlab("Time frame (years)") + ylab("Decay rate (yr-1)") +
  geom_point(color='green4')+
  geom_point(aes(100, mean(na.omit(k_table_Sg$k_100))), color='green4')+
  geom_errorbar(aes(100, ymin=mean(na.omit(k_table_Sg$k_100))-std(na.omit(k_table_Sg$k_100)), ymax=mean(na.omit(k_table_Sg$k_100))+std(na.omit(k_table_Sg$k_100))), color='green4')+
  geom_point(aes(100, mean(na.omit(k_table_Mg$k_100))), color='blue')+
  geom_errorbar(aes(100, ymin=mean(na.omit(k_table_Mg$k_100))-std(na.omit(k_table_Mg$k_100)), ymax=mean(na.omit(k_table_Mg$k_100))+std(na.omit(k_table_Mg$k_100))), color='blue')+
  geom_point(aes(100, mean(na.omit(k_table_Sm$k_100))), color='orange')+
  geom_errorbar(aes(100, ymin=mean(na.omit(k_table_Sm$k_100))-std(na.omit(k_table_Sm$k_100)), ymax=mean(na.omit(k_table_Sm$k_100))+std(na.omit(k_table_Sm$k_100))), color='orange')+
  
  geom_point(aes(150, mean(na.omit(k_table_Sg$k_150))), color='green4')+
  geom_errorbar(aes(150, ymin=mean(na.omit(k_table_Sg$k_150))-std(na.omit(k_table_Sg$k_150)), ymax=mean(na.omit(k_table_Sg$k_150))+std(na.omit(k_table_Sg$k_150))), color='green4')+
  geom_point(aes(150, mean(na.omit(k_table_Mg$k_150))), color='blue')+
  geom_errorbar(aes(150, ymin=mean(na.omit(k_table_Mg$k_150))-std(na.omit(k_table_Mg$k_150)), ymax=mean(na.omit(k_table_Mg$k_150))+std(na.omit(k_table_Mg$k_150))), color='blue')+  
  geom_point(aes(150, mean(na.omit(k_table_Sm$k_150))), color='orange')+
  geom_errorbar(aes(150, ymin=mean(na.omit(k_table_Sm$k_150))-std(na.omit(k_table_Sm$k_150)), ymax=mean(na.omit(k_table_Sm$k_150))+std(na.omit(k_table_Sm$k_150))), color='orange')+  
  
  geom_point(aes(300, mean(na.omit(k_table_Sg$k_300))), color='green4')+
  geom_errorbar(aes(300, ymin=mean(na.omit(k_table_Sg$k_300))-std(na.omit(k_table_Sg$k_300)), ymax=mean(na.omit(k_table_Sg$k_300))+std(na.omit(k_table_Sg$k_300))), color='green4')+
  geom_point(aes(300, mean(na.omit(k_table_Mg$k_300))), color='blue')+
  geom_errorbar(aes(300, ymin=mean(na.omit(k_table_Mg$k_300))-std(na.omit(k_table_Mg$k_300)), ymax=mean(na.omit(k_table_Mg$k_300))+std(na.omit(k_table_Mg$k_300))), color='blue')+
  geom_point(aes(300, mean(na.omit(k_table_Sm$k_300))), color='orange')+
  geom_errorbar(aes(300, ymin=mean(na.omit(k_table_Sm$k_300))-std(na.omit(k_table_Sm$k_300)), ymax=mean(na.omit(k_table_Sm$k_300))+std(na.omit(k_table_Sm$k_300))), color='orange')+
  
  geom_point(aes(500, mean(na.omit(k_table_Sg$k_500))), color='green4')+
  geom_errorbar(aes(500, ymin=mean(na.omit(k_table_Sg$k_500))-std(na.omit(k_table_Sg$k_500)), ymax=mean(na.omit(k_table_Sg$k_500))+std(na.omit(k_table_Sg$k_500))), color='green4')+
  geom_point(aes(500, mean(na.omit(k_table_Mg$k_500))), color='blue')+
  geom_errorbar(aes(500, ymin=mean(na.omit(k_table_Mg$k_500))-std(na.omit(k_table_Mg$k_500)), ymax=mean(na.omit(k_table_Mg$k_500))+std(na.omit(k_table_Mg$k_500))), color='blue')+
  geom_point(aes(500, mean(na.omit(k_table_Sm$k_500))), color='orange')+
  geom_errorbar(aes(500, ymin=mean(na.omit(k_table_Sm$k_500))-std(na.omit(k_table_Sm$k_500)), ymax=mean(na.omit(k_table_Sm$k_500))+std(na.omit(k_table_Sm$k_500))), color='orange')+
  
  geom_point(aes(1000, mean(na.omit(k_table_Sg$k_1000))), color='green4')+
  geom_errorbar(aes(1000, ymin=mean(na.omit(k_table_Sg$k_1000))-std(na.omit(k_table_Sg$k_1000)), ymax=mean(na.omit(k_table_Sg$k_1000))+std(na.omit(k_table_Sg$k_1000))), color='green4')+
  geom_point(aes(1000, mean(na.omit(k_table_Mg$k_1000))), color='blue')+
  geom_errorbar(aes(1000, ymin=mean(na.omit(k_table_Mg$k_1000))-std(na.omit(k_table_Mg$k_1000)), ymax=mean(na.omit(k_table_Mg$k_1000))+std(na.omit(k_table_Mg$k_1000))), color='blue')+
  geom_point(aes(1000, mean(na.omit(k_table_Sm$k_1000))), color='orange')+
  geom_errorbar(aes(1000, ymin=mean(na.omit(k_table_Sm$k_1000))-std(na.omit(k_table_Sm$k_1000)), ymax=mean(na.omit(k_table_Sm$k_1000))+std(na.omit(k_table_Sm$k_1000))), color='orange')+
  
  geom_point(aes(1500, mean(na.omit(k_table_Sg$k_1500))), color='green4')+
  geom_errorbar(aes(1500, ymin=mean(na.omit(k_table_Sg$k_1500))-std(na.omit(k_table_Sg$k_1500)), ymax=mean(na.omit(k_table_Sg$k_1500))+std(na.omit(k_table_Sg$k_1500))), color='green4')+
  geom_point(aes(1500, mean(na.omit(k_table_Mg$k_1500))), color='blue')+
  geom_errorbar(aes(1500, ymin=mean(na.omit(k_table_Mg$k_1500))-std(na.omit(k_table_Mg$k_1500)), ymax=mean(na.omit(k_table_Mg$k_1500))+std(na.omit(k_table_Mg$k_1500))), color='blue')+
  geom_point(aes(1500, mean(na.omit(k_table_Sm$k_1500))), color='orange')+
  geom_errorbar(aes(1500, ymin=mean(na.omit(k_table_Sm$k_1500))-std(na.omit(k_table_Sm$k_1500)), ymax=mean(na.omit(k_table_Sm$k_1500))+std(na.omit(k_table_Sm$k_1500))), color='orange')+
  
  geom_point(aes(2000, mean(na.omit(k_table_Sg$k_2000))), color='green4')+
  geom_errorbar(aes(2000, ymin=mean(na.omit(k_table_Sg$k_2000))-std(na.omit(k_table_Sg$k_2000)), ymax=mean(na.omit(k_table_Sg$k_2000))+std(na.omit(k_table_Sg$k_2000))), color='green4')+
  geom_point(aes(2000, mean(na.omit(k_table_Mg$k_2000))), color='blue')+
  geom_errorbar(aes(2000, ymin=mean(na.omit(k_table_Mg$k_2000))-std(na.omit(k_table_Mg$k_2000)), ymax=mean(na.omit(k_table_Mg$k_2000))+std(na.omit(k_table_Mg$k_2000))), color='blue')+
  geom_point(aes(2000, mean(na.omit(k_table_Sm$k_2000))), color='orange')+
  geom_errorbar(aes(2000, ymin=mean(na.omit(k_table_Sm$k_2000))-std(na.omit(k_table_Sm$k_2000)), ymax=mean(na.omit(k_table_Sm$k_2000))+std(na.omit(k_table_Sm$k_2000))), color='orange')+

  geom_point(data= k_rev,mapping = aes(k_rev$Max_Age, k_rev$k), color='green', shape= 17)




# results distribution

#get coordinates from B dataframe

mapa<-k_table[rowSums(is.na(k_table)) != ncol(k_table)-2,]

CoordR<-subset(B, B$Core %in% mapa$ID==TRUE)

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
       filename =  "Estimated k all map.jpg",
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

# model first 150 years Pb210 ---------------------------------------------------

TPbandC$FAge<-TPbandC$Age.Pb
Data_i<-subset(TPbandC, TPbandC$FAge < 150)
Data_t<-tendency(Data_i, pnames="150Pb")

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

fit_150Pb<-OCModel(DataAM, nwpath="Decay2023/150_Pb")

#eliminate some cores after visual check, we eliminate: Sg_316, Sg_317, Sg_321, Sm_004, Sm_022
fit_150Pb[c(26:28, 30,32), "k"]<-NA









# model first 100 years C14 ---------------------------------------------------

TPbandC$FAge<-TPbandC$Age.C
Data_i<-subset(TPbandC, TPbandC$FAge < 150)
Data_t<-tendency(Data_i, pnames="150_C")

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

fit_150_C<-OCModel(DataAM, nwpath="Decay2023/150_C")

#eliminate some cores after visual check, we eliminate: Sg_323
fit_150_C[c(22), "k"]<-NA


# comp 150 Pb vs C ---------------------------------------------------

PbvsC<-merge(fit_150Pb, fit_150_C, by = 'ID', all.x=T, all.y=F)

mPbvsC<-melt(PbvsC[,c(1,2,4,8)], id = c("ID","Ecosystem.x"))

mPbvsC$variable <- as.character(mPbvsC$variable)
mPbvsC$variable[mPbvsC$variable == 'k.x'] <- 'Modelo Pb'
mPbvsC$variable[mPbvsC$variable == 'k.y'] <- 'Modelo C'


ggplot(mPbvsC, aes(variable, value))+ ylab("Decay rate 150 yr (yr-1)")+
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem.x))

pairwise.wilcox.test(mPbvsC$value, mPbvsC$variable,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)



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
