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

B <- subset(Cores, Cores$V.vs.B != "Bare")
#Bare <- subset(A, A[,5]=='Bare')
length(unique(B$Core))
SingleCore<-B[!duplicated(B$Core),]

table(SingleCore$Ecosystem)


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

pairwise.wilcox.test(CM$Av_Mud, DT2$C_Gr,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

DT2<-DT
DT2$C_Gr <- recode(DT2$C_Gr, DEC = 'DEC',
                   INC  = 'N',
                   NT = 'N')

ggplot(CM, aes(DT2$C_Gr, CM$Av_Mud)) +
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



# chronological models estimation and data frame with age -----------------



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



TAll = filter(C, !is.na(Age) | !is.na(Age.Pb) | !is.na(Age.C))# df with core with any model
length(unique(TAll$Core))

TPb = subset(C, !is.na(Age) | !is.na(Age.Pb) )# df with core with Pb model
length(unique(TPb$Core))

TC = filter(C, !is.na(Age.C))# df with core with C model
length(unique(TC$Core))


# Correlation with time ---------------------------------------------------



#Only cores with Pb, if there is 14C with 14C, if not only Pb 


Folder = "Decay2023_Pb"
dir.create(Folder)

#### Homogenize Age


TPb$FAge<-"NA"
TPb$FAge<-as.numeric(TPb$FAge)

for (i in 1:nrow(TPb)) {
  
  if (is.na(TPb[i,"Age"]) == FALSE) {TPb[i,"FAge"]<-TPb[i,"Age"]} 
  
  else {
    
    if (is.na(TPb[i,"Age.Pb"]) == FALSE) {TPb[i,"FAge"]<-TPb[i,"Age.Pb"]}
    
    else {if (is.na(TPb[i,"Age.C"]) == FALSE) {TPb[i,"FAge"]<-TPb[i,"Age.C"]}}}}

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

ten_gr<-tendency(TPb, pnames="prueba")

TDEC<-ten_gr[[2]]

tendPb<-ten_gr[[1]]

#### Count cores per group and get the percentages

NGr <- tendPb %>% group_by(C_Gr) %>% count()
NGr %>% mutate(proc = ((n * 100) / sum(NGr[, 2])))

tendPb %>% group_by(Ecosystem, C_Gr) %>% count()
NGrSg <- subset(tendPb, Ecosystem == "Seagrass") %>% group_by(C_Gr) %>% count()
NGrSg %>% mutate(proc = ((n * 100) / sum(NGrSg[, 2])))

NGrSm <-
  subset(tendPb, Ecosystem == "Tidal Marsh") %>% group_by(C_Gr) %>% count()
NGrSm %>% mutate(proc = ((n * 100) / sum(NGrSm[, 2])))

NGrMg <- subset(tendPb, Ecosystem == "Mangrove") %>% group_by(C_Gr) %>% count()
NGrMg %>% mutate(proc = ((n * 100) / sum(NGrMg[, 2])))



# NLM time-Accumulated Mass -----------------------------------------------


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

max_depth<-OCModel(DataAM, nwpath="Decay2023_Pb/Max_Depth")

#eliminate outlayers:
#eliminate empty cores (no model): 
#eliminate some cores after visual check, we eliminate: Sg_041, Sg_310, Sg_316, Sg_321, Sm_010, Sm_022, Sm_68, Sm_69, 
#Sm_092, Sm_097, Sm_105
max_depth <- max_depth[-c(14, 18, 37, 41, 43, 46, 47, 60, 61, 66, 68, 69), ]



ggplot(max_depth, aes(x = k)) +
  geom_histogram()

shapiro.test(max_depth$k) #(>0.05 normal, <0.05 no normal)

# check if time frame has an effect 

ggplot(max_depth, aes(Max.Age, k))+
  geom_point()



# model first 100 years ---------------------------------------------------

Data_i<-subset(TPb, TPb$FAge < 100)
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

fit_100Pb<-OCModel(DataAM, nwpath="Decay2023_Pb/100")

#eliminate some cores after visual check, we eliminate: Sg_081, Sg_111, Sg_192, Sg_316, Sg_321, Sg_323, Sg_332, Sm_004, 
#Sm_068, Sm_069, Sm_092, Sm_097, Sm_105
fit_100Pb[c(10, 13, 21, 31, 33, 34, 35, 37, 49, 50, 55, 57, 58), "k"]<-NA


ggplot(fit_100Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

ggplot(fit_100Pb, aes( Max.Age, k))+
  geom_point()


pairwise.wilcox.test(fit_100Pb$k, fit_100Pb$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)



# model 100-150 years ---------------------------------------------------

Data_i<-subset(TPb, TPb$FAge < 150)
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

fit_150Pb<-OCModel(DataAM, MA= 100, nwpath="Decay2023_Pb/150")


#eliminate some cores after visual check, we eliminate: Mg_023, Sm_004
fit_150Pb <- fit_150Pb[-c(6,40), ]

#eliminate empty cores (no model)
fit_150Pb<-fit_150Pb[!is.na(fit_150Pb$P),]


ggplot(fit_150Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(fit_150Pb$k, fit_150Pb$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model 150-300 years ---------------------------------------------------


Data_i<-subset(TPb, TPb$FAge < 300)
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

fit_300Pb<-OCModel(DataAM, MA= 150, nwpath="Decay2023_Pb/300")


#eliminate outlayers: 

ggplot(fit_300Pb, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: Sg_179, Sg_311, Sm_004
fit_300Pb <- fit_300Pb[-c( 28, 37, 45), ]

#eliminate empty cores (no model)
fit_300Pb<-fit_300Pb[!is.na(fit_300Pb$P),]


ggplot(fit_300Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()


pairwise.wilcox.test(fit_300Pb$k, fit_300$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


# model 300-500 years ---------------------------------------------------

Data_i<-subset(TPb, TPb$FAge < 500)
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

fit_500Pb<-OCModel(DataAM, MA= 300, nwpath="Decay2023_Pb/500")


#eliminate some cores after visual check, we eliminate: Sg_019, Sg_041, Sg_190, Sg_192, Sg_310, Sg_312, Sg_314, Sg_323, Sm_004
fit_500Pb <- fit_500Pb[-c(12, 14, 29, 31, 37, 39, 40, 44, 46 ), ]

#eliminate empty cores (no model)
fit_50Pb0<-fit_500Pb[!is.na(fit_500Pb$P),]


ggplot(fit_500Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(fit_500Pb$k, fit_500Pb$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model 500-1000 years ---------------------------------------------------


Data_i<-subset(TPb, TPb$FAge < 1000)
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

fit_1000Pb<-OCModel(DataAM, MA= 500, nwpath="Decay2023_Pb/1000")

#eliminate some cores after visual check, we eliminate: Sg_019, Sg_041, Sg_111, Sg_117, Sg_195, Sg_312, Sg_314, Sm_004, Sm_010
fit_1000Pb <- fit_1000Pb[-c( 12, 13, 19, 21, 32, 40, 41, 47, 48), ]

#eliminate empty cores (no model)
fit_1000Pb<-fit_1000Pb[!is.na(fit_1000Pb$P),]


ggplot(fit_1000Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()


pairwise.wilcox.test(fit_1000Pb$k, fit_1000Pb$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model 1000-1500 years ---------------------------------------------------


Data_i<-subset(TPb, TPb$FAge < 1500)
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

fit_1500Pb<-OCModel(DataAM, MA= 1000, nwpath="Decay2023_Pb/1500")


#eliminate outlayers: 

ggplot(fit_1500Pb, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate: 
fit_1500Pb <- fit_1500Pb[-c(13, 17, 19, 20, 23, 25, 31, 39, 40, 46, 48 ), ]

#eliminate empty cores (no model)
fit_1500Pb<-fit_1500Pb[!is.na(fit_1500Pb$P),]


ggplot(fit_1500Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()

pairwise.wilcox.test(fit_1500Pb$k, fit_1500Pb$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


# model 1500-2000 years ---------------------------------------------------


Data_i<-subset(TPb, TPb$FAge < 2000)
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

fit_2000Pb<-OCModel(DataAM, MA= 1500, nwpath="Decay2023_Pb/2000")


#eliminate outlayers: 

ggplot(fit_2000Pb, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate:
fit_2000Pb <- fit_2000Pb[-c( 29), ]

#eliminate empty cores (no model)
fit_2000Pb<-fit_2000Pb[!is.na(fit_2000Pb$P),]


ggplot(fit_2000Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()


pairwise.wilcox.test(fit_2000Pb$k, fit_2000Pb$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

# model > 2000 years ---------------------------------------------------

Data_i<-TPb
Data_t<-tendency(Data_i, pnames=">2000_Pb")

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

fit_m2000Pb<-OCModel(DataAM, MA= 2000, nwpath="Decay2023_Pb/more_2000")


#eliminate outlayers: 

ggplot(fit_m2000Pb, aes(x = k)) +
  geom_histogram()

#eliminate some cores after visual check, we eliminate:
fit_m2000Pb <- fit_m2000Pb[-c(18, 40), ]

#eliminate empty cores (no model)
fit_m2000Pb<-fit_m2000Pb[!is.na(fit_m2000Pb$P),]


ggplot(fit_m2000Pb, aes( Ecosystem, k))+
  geom_boxplot()+
  geom_jitter()







# Final table and plots -------------------------------------------------------

k_tablePb <-merge(fit_100Pb[,c(1,4)], fit_150Pb[,c(1,4)], by = "ID", all = TRUE)
k_tablePb <-merge(k_tablePb, fit_300Pb[,c(1,4)], by = "ID", all = TRUE)
k_tablePb <-merge(k_tablePb, fit_500Pb[,c(1,4)], by = "ID", all = TRUE)
k_tablePb <-merge(k_tablePb, fit_1000Pb[,c(1,4)], by = "ID", all = TRUE)
k_tablePb <-merge(k_tablePb, fit_1500Pb[,c(1,4)], by = "ID", all = TRUE)
k_tablePb <-merge(k_tablePb, fit_2000Pb[,c(1,4)], by = "ID", all = TRUE)
k_tablePb <-merge(k_tablePb, fit_m2000Pb[,c(1,4, 5)], by = "ID", all = TRUE)

colnames(k_tablePb)<-c("ID", "k_100", "k_150", "k_300", "k_500", "k_1000","k_1500", "k_2000", "k_m2000", "Max_Age")
k_tablePb$k_100<-as.numeric(k_tablePb$k_100)
k_tablePb$k_300<-as.numeric(k_tablePb$k_300)

names(SingleCore)[names(SingleCore) == 'Core'] <- 'ID'

k_tablePb<-merge(k_tablePb, SingleCore[,c(1, 3, 7, 6, 8, 11, 12)], by = 'ID', all.x=T, all.y=F)
k_tablePb<-merge(k_tablePb, SAR[,c(1,5)], by = 'ID', all.x=T, all.y=F)
k_tablePb<-merge(k_tablePb, CM[,c(1,9,11)], by = 'ID', all.x=T, all.y=F)



#normal distribution and significant differences among ecosystems

shapiro.test(k_tablePb$k_100) #normal if pvalue > than 0.05

apply(k_tablePb[,c(2:8)], FUN=shapiro.test, MARGIN = 2)


# differences among species

pairwise.wilcox.test(k_tablePb$k_100, k_tablePb$Specie,
                     p.adjust.method = "BH")




# correlation decay rate (150) and SAR ------------------------------------


plot(k_tablePb_Sg$k_100, k_tablePb_Sg$Av_Mud_25)

ggplot(k_tablePb,aes(k_1000, SAR))+
  geom_point(aes(color=Bioregions))

cor.test(k_tablePb$k_150, k_tablePb$SAR, method=c("pearson"))

ggplot(k_tablePb, aes(k_150, SAR))+
  geom_point(aes(color=Ecosystem))

shapiro.test(k_tablePb$k_150)
shapiro.test(k_tablePb$Av_Mud_25)


# correlation decay rate and Mud ------------------------------------

# first we estimate the average mud content for the studied time frame
# we use the max age of each core

x<-split(TPb, TPb$Core)

Mud150_a<-fit_150Pb
Mud150_a$Mud<-"NA"

for (i in 1:nrow(Mud150_a)) {
  
  data<-as.data.frame(x[Mud150_a[i,"ID"]])
  colnames(data)<-colnames(TPb)
  data_a<-subset(data, data$FAge<=Mud150_a[i, "Max.Age"])
  Mud150_a[i, "Mud"]<-mean(data_a$Mud, na.rm=TRUE)}

Mud150_a$Mud<-as.numeric(Mud150_a$Mud)

Mud1000_a<-fit_1000Pb
Mud1000_a$Mud<-"NA"

for (i in 1:nrow(Mud1000_a)) {
  
  data<-as.data.frame(x[Mud1000_a[i,"ID"]])
  colnames(data)<-colnames(TPb)
  data_a<-subset(data, data$FAge<=Mud1000_a[i, "Max.Age"])
  Mud1000_a[i, "Mud"]<-mean(data_a$Mud, na.rm=TRUE)}

Mud1000_a$Mud<-as.numeric(Mud1000_a$Mud)

shapiro.test(Mud150_a$k)
shapiro.test(as.numeric(Mud150_a$Mud))
plot(Mud150_a$k, Mud150_a$Mud)
cor.test(as.numeric(Mud150_a$k), as.numeric(Mud150_a$Mud), method=c("spearman"))

ggplot(Mud150_a, aes(k, Mud))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c('blue', 'green4', "orange"))


plot(Mud1000_a$k, Mud1000_a$Mud)
cor.test(as.numeric(Mud1000_a$k), as.numeric(Mud1000_a$Mud), method=c("spearman"))

plot(k_tablePb_Sg$k_150, k_tablePb_Sg$Av_Mud_25)

ggplot(k_tablePb, aes(k_150, Av_Mud_25))+
  geom_point(aes(color=Ecosystem))


# Mud distribution of fitted cores vs mud distribution global

library(moments)#skewness


# whole core collection
histogram(CM$Av_Mud_25, breaks = 10)
mean(CM$Av_Mud_25, na.rm=TRUE)
std(CM$Av_Mud_25)
nrow(CM[complete.cases(CM[,"Av_Mud_25"]),])
skewness(CM$Av_Mud_25, na.rm = TRUE)

P1<-ggplot(CM, aes(x=Av_Mud_25))+ ylab("Whole core collection")+
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 10)+
  scale_y_continuous(labels = scales::percent)+
  annotate("text", x=60, y=0.18, label= "27.2±0.95 (mean±SE)")+
  annotate("text", x=60, y=0.15, label= "skewness=1.1")+
  annotate("text", x=60, y=0.12, label= "n=248")+
  theme(axis.title.x=element_blank())

# cores fitted to last 100-150 yr model
Mud150<-k_tablePb[complete.cases(k_tablePb[,"k_150"]),]

histogram(Mud150$Av_Mud_25, breaks = 10)
mean(Mud150$Av_Mud_25, na.rm=TRUE)
std(Mud150$Av_Mud_25)
nrow(Mud150[complete.cases(Mud150[,"Av_Mud_25"]),])
skewness(Mud150$Av_Mud_25, na.rm = TRUE)

P2<-ggplot(Mud150, aes(x=Av_Mud_25))+ ylab("Cores 100-150 yr model") + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 10)+
  scale_y_continuous(labels = scales::percent)+
  annotate("text", x=60, y=0.25, label= "23.5±4.8 (mean±SE)")+
  annotate("text", x=60, y=0.2, label= "skewness=1.3")+
  annotate("text", x=60, y=0.16, label= "n=21")+
  theme(axis.title.x=element_blank())

# cores fitted to last 500-1000 yr model
Mud1000<-k_tablePb[complete.cases(k_tablePb[,"k_1000"]),]

histogram(Mud1000$Av_Mud_25, breaks = 10)
mean(Mud1000$Av_Mud_25, na.rm=TRUE)
std(Mud1000$Av_Mud_25)
nrow(Mud1000[complete.cases(Mud1000[,"Av_Mud_25"]),])
skewness(Mud1000$Av_Mud_25, na.rm = TRUE)

P3<-ggplot(Mud1000, aes(x=Av_Mud_25))+ ylab("Cores 500-1000 yr model") + xlab("Mud (<0.063 mm) concentration (%)") +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 10)+
  scale_y_continuous(labels = scales::percent)+
  annotate("text", x=60, y=0.2, label= "30.9±6.4 (mean±SE)")+
  annotate("text", x=60, y=0.17, label= "skewness=1")+
  annotate("text", x=60, y=0.14, label= "n=22")


gridExtra::grid.arrange(P1, P2, P3, ncol = 1)

### boxplots and wilcox

CM$label<-"All"
Mud150$label<-"k 100-150"
Mud1000$label<-"k 500-1000"


ggplot(CM, aes(label,Av_Mud_25))+ ylab("Mud % (<0.063 mm)") +
  geom_boxplot()+
  geom_boxplot(data=Mud150, aes(label,Av_Mud_25))+
  geom_boxplot(data=Mud1000, aes(label,Av_Mud_25))+
  geom_jitter(aes(color=Ecosystem))+
  geom_jitter(data=Mud150, aes(color=Ecosystem))+
  geom_jitter(data=Mud1000, aes(color=Ecosystem))+
  scale_color_manual(values=c('blue', 'green4', "orange"))
  

temp<-rbind(CM[,c(1,11,14)], Mud150[,c(1,19,20)], Mud1000[,c(1,19,20)])
temp$Av_Mud_25<-as.numeric(temp$Av_Mud_25)

pairwise.wilcox.test(temp$Av_Mud_25, temp$label,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)


global_mud<-ggplot(CM, aes(Ecosystem,Av_Mud_25))+ ylab("Mud % top 25 cm (<0.063 mm)") +
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))+
  scale_color_manual(values=c('blue', 'green4', "orange"))+
  ylim(0, 110)


pairwise.wilcox.test(CM$Av_Mud_25, CM$Ecosystem,
                     p.adjust.method = "BH") # are significantly different (p < 0.05)

ggsave(
  plot = global_mud,
  path = Folder,
  filename = "global_mud.jpg",
  units = "cm",
  width = 12,
  height = 7
)


#only for seagrass

CM_Sg<-subset(CM, Ecosystem == "Seagrass")

mean(CM_Sg$Av_Mud_25, na.rm=TRUE)
std(CM_Sg$Av_Mud_25)


#check for outlayers (https://www.r-bloggers.com/2016/12/outlier-detection-and-treatment-with-r/)

# visual check
boxplot(k_100 ~ Ecosystem, data=k_tablePb, main="Decay 150yr by Ecosystem")  # clear pattern is noticeable.
boxplot(SAR ~ Ecosystem, data=k_tablePb, main="SAR by Ecosystem")  # this may not be significant, as day of week variable is a subset of the month var.

#Cook’s Distance
mod <- lm(SAR ~ k_100, data=k_tablePb)
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

# one outlayer: 100

cor.test(k_tablePb[,-c(51,100)]$k_100, k_tablePb[,-c(51,100)]$SAR, method=c("pearson"))

ggplot(k_tablePb[,-c(51,100)],aes(k_100, SAR))+
  geom_point(aes(color=Ecosystem))








# summary table (manuscript Table 1) -----------------------------------


k_tablePb_Mg<-subset(k_tablePb, Ecosystem=='Mangrove')
k_tablePb_Mg<-sapply(k_tablePb_Mg,FUN=as.numeric)
k_tablePb_Sg<-subset(k_tablePb, Ecosystem=='Seagrass')
k_tablePb_Sg<-sapply(k_tablePb_Sg,FUN=as.numeric)
k_tablePb_Sm<-subset(k_tablePb, Ecosystem=='Tidal Marsh')
k_tablePb_Sm<-sapply(k_tablePb_Sm,FUN=as.numeric)

std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
cnt <- function(x) sum(!is.na(x))   

sum_table<-as.data.frame(colMeans(k_tablePb_Mg[,c(2:8)], na.rm = TRUE))
sum_table[,2]<-as.data.frame(apply(k_tablePb_Mg[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,3]<-as.data.frame(apply(k_tablePb_Mg[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,4]<-colMeans(k_tablePb_Sg[,c(2:8)], na.rm = TRUE)
sum_table[,5]<-as.data.frame(apply(k_tablePb_Sg[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,6]<-as.data.frame(apply(k_tablePb_Sg[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,7]<-colMeans(k_tablePb_Sm[,c(2:8)], na.rm = TRUE)
sum_table[,8]<-as.data.frame(apply(k_tablePb_Sm[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,9]<-as.data.frame(apply(k_tablePb_Sm[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,10]<-colMeans(k_tablePb[,c(2:8)], na.rm = TRUE)
sum_table[,11]<-as.data.frame(apply(k_tablePb[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,12]<-as.data.frame(apply(k_tablePb[,c(2:8)], FUN=cnt, MARGIN = 2))

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
k_revS<- k_rev [-c(1,6,8,9),]
k_revM<- k_rev [c(8,9),]



# boxplot by timeframe figure

mk_tablePb<-melt(k_tablePb[,c(1:8, 11, 12)], id = c("ID","Ecosystem", "Specie")) 

mk_tablePb$variable <- as.character(mk_tablePb$variable)
mk_tablePb$variable[mk_tablePb$variable == 'k_100'] <- '0-100 yr'
mk_tablePb$variable[mk_tablePb$variable == 'k_150'] <- '100-150 yr'
mk_tablePb$variable[mk_tablePb$variable == 'k_300'] <- '150-300 yr'
mk_tablePb$variable[mk_tablePb$variable == 'k_500'] <- '300-500 yr'
mk_tablePb$variable[mk_tablePb$variable == 'k_1000'] <- '500-1000 yr'
mk_tablePb$variable[mk_tablePb$variable == 'k_1500'] <- '1000-1500 yr'
mk_tablePb$variable[mk_tablePb$variable == 'k_2000'] <- '1500-2000 yr'

mk_tablePb$value<-as.numeric(mk_tablePb$value)

ggplot(transform(mk_tablePb,
                 variable=factor(variable,levels=c('0-100 yr','100-150 yr','150-300 yr', '300-500 yr', '500-1000 yr', '1000-1500 yr', '1500-2000 yr'))),
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


mk_tablePb2<-subset(mk_tablePb, mk_tablePb$variable == "100-150 yr" | mk_tablePb$variable == "500-1000 yr")

box_150_1000<-  ggplot(transform(mk_tablePb2,
                 variable=factor(variable,levels=c('100-150 yr',"500-1000 yr"))),
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
  plot = box_150_1000,
  path = Folder,
  filename = "box_150_1000.jpg",
  units = "cm",
  width = 12,
  height = 7
)


ggplot(transform(mk_tablePb2,
                 variable=factor(variable,levels=c('100-150 yr',"500-1000 yr"))),
       aes(Ecosystem, value))+ ggtitle("Decay rates by ecosystem and time frame")+ ylab("Decay rate (yr-1)") +
  geom_boxplot()+
  geom_jitter(aes(color=Specie))+
  facet_wrap(~variable)+
  #scale_color_manual(values=c('blue', 'green4', "orange"))+
  theme(#legend.position = c(1, 0),
    #legend.justification = c(1, 0), 
    axis.title.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),)






# exponential model to predict k  -----------------------------------------

#fitting table

f_table_Sg<-as.data.frame(colMeans(k_tablePb[c(7:51),c(2:8)], na.rm = TRUE))
f_table_Sg[,2]<-c(50, 125, 225, 400, 750, 1250, 1750)
f_table_Sg[c(8:11),1]<-na.omit(k_tablePb[,9])
f_table_Sg[c(8:11),2]<-na.omit(k_tablePb[,10])
f_table_Sg[c(12:16), c(1,2)]<-k_revS[,c(2:3)]
f_table[c(17:18), c(1,2)]<-k_revM[,c(2:3)]

f_table_Sg<-f_table[-9,]
colnames(f_table_Sg)<-c("k","timeframe")

f_table_Tm<-as.data.frame(colMeans(k_tablePb[c(52:76),c(2:7)], na.rm = TRUE))
f_table_Tm[,2]<-c(50, 125, 225, 400, 750, 1250)
f_table_Tm[c(7:8), c(1,2)]<-k_revM[,c(2:3)]
colnames(f_table_Tm)<-c("k","timeframe")

f_table_Mg<-as.data.frame(colMeans(k_tablePb[c(1:6),c(2:6)], na.rm = TRUE))
f_table_Mg[,2]<-c(50, 125, 225, 400, 750)
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
fitSg[, 2] <- kchange(c(1:5000), 0.0237, -0.002)
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
fitTm[, 2] <- kchange(c(1:5000), 0.021, -0.003)
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
fitMg[, 2] <- kchange(c(1:5000), 0.0268, -0.0038)
colnames(fitMg) <- list("timeframe", "predict")



# k vs time frame fitting figure ------------------------------------------


std <- function(x) sd(x)/sqrt(length(x))

fit_figPb<-
  ggplot(k_tablePb_Sg, aes( Max_Age, k_m2000))+ ggtitle("Decay rate by time frame") + xlab("Time frame (years)") + ylab("Decay rate (yr-1)") +
  geom_point(color='green4')+
  geom_point(aes(50, mean(na.omit(k_tablePb_Sg$k_100))), color='green4')+
  geom_errorbar(aes(50, ymin=mean(na.omit(k_tablePb_Sg$k_100))-std(na.omit(k_tablePb_Sg$k_100)), ymax=mean(na.omit(k_tablePb_Sg$k_100))+std(na.omit(k_tablePb_Sg$k_100))), color='green4')+
  geom_point(aes(50, mean(na.omit(k_tablePb_Mg$k_100))), color='blue')+
  geom_errorbar(aes(50, ymin=mean(na.omit(k_tablePb_Mg$k_100))-std(na.omit(k_tablePb_Mg$k_100)), ymax=mean(na.omit(k_tablePb_Mg$k_100))+std(na.omit(k_tablePb_Mg$k_100))), color='blue')+
  geom_point(aes(50, mean(na.omit(k_tablePb_Sm$k_100))), color='orange')+
  geom_errorbar(aes(50, ymin=mean(na.omit(k_tablePb_Sm$k_100))-std(na.omit(k_tablePb_Sm$k_100)), ymax=mean(na.omit(k_tablePb_Sm$k_100))+std(na.omit(k_tablePb_Sm$k_100))), color='orange')+
  
  geom_point(aes(125, mean(na.omit(k_tablePb_Sg$k_150))), color='green4')+
  geom_errorbar(aes(125, ymin=mean(na.omit(k_tablePb_Sg$k_150))-std(na.omit(k_tablePb_Sg$k_150)), ymax=mean(na.omit(k_tablePb_Sg$k_150))+std(na.omit(k_tablePb_Sg$k_150))), color='green4')+
  geom_point(aes(125, mean(na.omit(k_tablePb_Mg$k_150))), color='blue')+
  geom_errorbar(aes(125, ymin=mean(na.omit(k_tablePb_Mg$k_150))-std(na.omit(k_tablePb_Mg$k_150)), ymax=mean(na.omit(k_tablePb_Mg$k_150))+std(na.omit(k_tablePb_Mg$k_150))), color='blue')+  
  geom_point(aes(125, mean(na.omit(k_tablePb_Sm$k_150))), color='orange')+
  geom_errorbar(aes(125, ymin=mean(na.omit(k_tablePb_Sm$k_150))-std(na.omit(k_tablePb_Sm$k_150)), ymax=mean(na.omit(k_tablePb_Sm$k_150))+std(na.omit(k_tablePb_Sm$k_150))), color='orange')+  
  
  geom_point(aes(175, mean(na.omit(k_tablePb_Sg$k_300))), color='green4')+
  geom_errorbar(aes(175, ymin=mean(na.omit(k_tablePb_Sg$k_300))-std(na.omit(k_tablePb_Sg$k_300)), ymax=mean(na.omit(k_tablePb_Sg$k_300))+std(na.omit(k_tablePb_Sg$k_300))), color='green4')+
  geom_point(aes(175, mean(na.omit(k_tablePb_Mg$k_300))), color='blue')+
  geom_errorbar(aes(175, ymin=mean(na.omit(k_tablePb_Mg$k_300))-std(na.omit(k_tablePb_Mg$k_300)), ymax=mean(na.omit(k_tablePb_Mg$k_300))+std(na.omit(k_tablePb_Mg$k_300))), color='blue')+
  geom_point(aes(175, mean(na.omit(k_tablePb_Sm$k_300))), color='orange')+
  geom_errorbar(aes(175, ymin=mean(na.omit(k_tablePb_Sm$k_300))-std(na.omit(k_tablePb_Sm$k_300)), ymax=mean(na.omit(k_tablePb_Sm$k_300))+std(na.omit(k_tablePb_Sm$k_300))), color='orange')+
  
  geom_point(aes(400, mean(na.omit(k_tablePb_Sg$k_500))), color='green4')+
  geom_errorbar(aes(400, ymin=mean(na.omit(k_tablePb_Sg$k_500))-std(na.omit(k_tablePb_Sg$k_500)), ymax=mean(na.omit(k_tablePb_Sg$k_500))+std(na.omit(k_tablePb_Sg$k_500))), color='green4')+
  geom_point(aes(400, mean(na.omit(k_tablePb_Mg$k_500))), color='blue')+
  geom_errorbar(aes(400, ymin=mean(na.omit(k_tablePb_Mg$k_500))-std(na.omit(k_tablePb_Mg$k_500)), ymax=mean(na.omit(k_tablePb_Mg$k_500))+std(na.omit(k_tablePb_Mg$k_500))), color='blue')+
  geom_point(aes(400, mean(na.omit(k_tablePb_Sm$k_500))), color='orange')+
  geom_errorbar(aes(400, ymin=mean(na.omit(k_tablePb_Sm$k_500))-std(na.omit(k_tablePb_Sm$k_500)), ymax=mean(na.omit(k_tablePb_Sm$k_500))+std(na.omit(k_tablePb_Sm$k_500))), color='orange')+
  
  geom_point(aes(750, mean(na.omit(k_tablePb_Sg$k_1000))), color='green4')+
  geom_errorbar(aes(750, ymin=mean(na.omit(k_tablePb_Sg$k_1000))-std(na.omit(k_tablePb_Sg$k_1000)), ymax=mean(na.omit(k_tablePb_Sg$k_1000))+std(na.omit(k_tablePb_Sg$k_1000))), color='green4')+
  geom_point(aes(750, mean(na.omit(k_tablePb_Mg$k_1000))), color='blue')+
  geom_errorbar(aes(750, ymin=mean(na.omit(k_tablePb_Mg$k_1000))-std(na.omit(k_tablePb_Mg$k_1000)), ymax=mean(na.omit(k_tablePb_Mg$k_1000))+std(na.omit(k_tablePb_Mg$k_1000))), color='blue')+
  geom_point(aes(750, mean(na.omit(k_tablePb_Sm$k_1000))), color='orange')+
  geom_errorbar(aes(750, ymin=mean(na.omit(k_tablePb_Sm$k_1000))-std(na.omit(k_tablePb_Sm$k_1000)), ymax=mean(na.omit(k_tablePb_Sm$k_1000))+std(na.omit(k_tablePb_Sm$k_1000))), color='orange')+
  
  geom_point(aes(1250, mean(na.omit(k_tablePb_Sg$k_1500))), color='green4')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(k_tablePb_Sg$k_1500))-std(na.omit(k_tablePb_Sg$k_1500)), ymax=mean(na.omit(k_tablePb_Sg$k_1500))+std(na.omit(k_tablePb_Sg$k_1500))), color='green4')+
  geom_point(aes(1250, mean(na.omit(k_tablePb_Mg$k_1500))), color='blue')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(k_tablePb_Mg$k_1500))-std(na.omit(k_tablePb_Mg$k_1500)), ymax=mean(na.omit(k_tablePb_Mg$k_1500))+std(na.omit(k_tablePb_Mg$k_1500))), color='blue')+
  geom_point(aes(1250, mean(na.omit(k_tablePb_Sm$k_1500))), color='orange')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(k_tablePb_Sm$k_1500))-std(na.omit(k_tablePb_Sm$k_1500)), ymax=mean(na.omit(k_tablePb_Sm$k_1500))+std(na.omit(k_tablePb_Sm$k_1500))), color='orange')+
  
  geom_point(aes(1750, mean(na.omit(k_tablePb_Sg$k_2000))), color='green4')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(k_tablePb_Sg$k_2000))-std(na.omit(k_tablePb_Sg$k_2000)), ymax=mean(na.omit(k_tablePb_Sg$k_2000))+std(na.omit(k_tablePb_Sg$k_2000))), color='green4')+
  geom_point(aes(1750, mean(na.omit(k_tablePb_Mg$k_2000))), color='blue')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(k_tablePb_Mg$k_2000))-std(na.omit(k_tablePb_Mg$k_2000)), ymax=mean(na.omit(k_tablePb_Mg$k_2000))+std(na.omit(k_tablePb_Mg$k_2000))), color='blue')+
  geom_point(aes(1750, mean(na.omit(k_tablePb_Sm$k_2000))), color='orange')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(k_tablePb_Sm$k_2000))-std(na.omit(k_tablePb_Sm$k_2000)), ymax=mean(na.omit(k_tablePb_Sm$k_2000))+std(na.omit(k_tablePb_Sm$k_2000))), color='orange')+
  
  geom_point(data= k_revS,mapping = aes(k_revS$Max_Age, k_revS$k), color='green', shape= 17)+
  geom_point(data= k_revM,mapping = aes(k_revM$Max_Age, k_revM$k), color='orange4', shape= 17)+
  
  geom_line(data= fitSg, mapping = aes(timeframe, predict), color='green4')+
  geom_line(data= fitTm, mapping = aes(timeframe, predict), color='orange')+
  geom_line(data= fitMg, mapping = aes(timeframe, predict), color='blue')+
  
  xlim(0,4000)+
  
    annotate("text", x=1500, y=0.025, color= "blue",  size = 5, label= expression(y == 0.027 * e ** (-0.004 * 
                                                                                             x)))+
    annotate("text", x=1500, y=0.02, color= "green4",size = 5,label= expression(y == 0.024 * e ** (-0.002 * 
                                                                              x)))+
  annotate("text", x=1500, y=0.015, color= "orange",size = 5,label= expression(y == 0.021 * e ** (-0.003 * 
                                                                          x)))
  

ggsave(
  plot = fit_figPb,
  path = Folder,
  filename = "fit_plot.jpg",
  units = "cm",
  width = 12,
  height = 7
)



# # results distribution (Fig 4) --------------------------------------------------


#get coordinates from B dataframe

WM <- map_data("world")

mapa<-k_tablePb[rowSums(is.na(k_tablePb)) != ncol(k_tablePb)-2,]

CoordR<-subset(B, B$Core %in% mapa$ID==TRUE)
CoordR$Lat<-as.numeric(CoordR$Lat)
CoordR$Long<-as.numeric(CoordR$Long)

CoordR %>%
  ggplot() + ggtitle("Estimated k distribution") + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
  #geom_point(aes(x = Long, y = Lat))+
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



# age and decay rates at 100 cm -----------------------------------------------------------


extract_depth_by_age<- function (df, AGE = 100) {
  
  core <- df[1,"Core"]
  c_length <- max(df[,"Max.Depth"])
  age_depth <- df[which.min(abs(AGE-df$Max.Depth)),"FAge"]
  
  output <- data.frame(core = core, Length= c_length, age_depth = age_depth)
  
  return(output)
  
}

x<-split(TPb, TPb$Core)

core_100_age_l <- lapply(X = x,  extract_depth_by_age, AGE = 100) # return a list
core_100_age <- as.data.frame(do.call(rbind, core_100_age_l)) # from list to dataframe

core_100_age <- subset(core_100_age, core_100_age$Length > 90)

mean(core_100_age$age_depth)
std(core_100_age$age_depth)


# from time to decay rate

#seagrass
0.024 * exp(-0.002 * 1953.206)

0.024 * exp(-0.002 * 1504.5)
0.024 * exp(-0.002 * 2401.9)

#mangrove
0.027 * exp(-0.004 * 1953.206)

0.027 * exp(-0.004 * 1504.5)
0.027 * exp(-0.004 * 2401.9)

#tidal marshes
0.021 * exp(-0.003 * 1953.206)

0.021 * exp(-0.003 * 1504.5)
0.021 * exp(-0.003 * 2401.9)
