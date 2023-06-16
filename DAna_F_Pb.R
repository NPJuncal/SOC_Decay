

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


pairwise.wilcox.test(fit_1000$k, fit_1000$Ecosystem,
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



#normal distribution and significant differences among ecosystems

shapiro.test(k_tablePb$k_100) #normal if pvalue > than 0.05

apply(k_tablePb[,c(2:8)], FUN=shapiro.test, MARGIN = 2)


# summary table (manuscript Table 1)

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
k_rev<- k_rev [-c(1,6),]

# boxplot by timeframe figure

mk_tablePb<-melt(k_tablePb[,c(1:8, 11)], id = c("ID","Ecosystem")) 

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

ggplot(transform(mk_tablePb2,
                 variable=factor(variable,levels=c('100-150 yr',"500-1000 yr"))),
       aes(Ecosystem, value))+ ggtitle("Decay rates by ecosystem and time frame")+ ylab("Decay rate (yr-1)") +
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))+
  facet_wrap(~variable)+
  scale_color_manual(values=c('blue', 'green4', "orange"))+
  theme(#legend.position = c(1, 0),
        #legend.justification = c(1, 0), 
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),)


#fitting table

f_table<-as.data.frame(colMeans(k_tablePb[,c(2:8)], na.rm = TRUE))
f_table[,2]<-c(50, 125, 225, 400, 750, 1250, 1750)
f_table[c(8:11),1]<-na.omit(k_tablePb[,9])
f_table[c(8:11),2]<-na.omit(k_tablePb[,10])
f_table[c(12:15), c(1,2)]<-k_rev[,c(2:3)]

f_table<-f_table[-9,]

# fit function k-timeframe

### exponential model to predict k in Posidonia meadows

colnames(f_table)<-c("k","timeframe")

kchange <- function(Tframe, A, C)
  (A * exp(C * Tframe))

model1 <-
  nls(
    k ~ kchange(timeframe, myA, myC),
    data = f_table,
    start = list(myA = 0.02, myC = 0.0003)
  )

# model 2 without 1500-2000 fittings
model2 <-
  nls(
    k ~ kchange(timeframe, myA, myC),
    data = f_table[-7,],
    start = list(myA = 0.02, myC = 0.0003)
  )



summary(model1)
plot(nlsResiduals(model1))

fitY1 <- kchange(c(1:2100), 0.027, -0.0018)

fitY1 <- as.data.frame(c(1:5000))
fitY1['new_col'] <- NA
fitY1[, 2] <- kchange(c(1:5000), 0.027, -0.0018)
colnames(fitY1) <- list("timeframe", "predict")


fitY2 <- as.data.frame(c(1:5000))
fitY2['new_col'] <- NA
fitY2[, 2] <- kchange(c(1:5000), 0.028, -0.0019)
colnames(fitY2) <- list("timeframe", "predict")


# correlation decay rate (150) and SAR ------------------------------------


plot(k_tablePb_Sg$k_100, k_tablePb_Sg$SAR)

ggplot(k_tablePb_Sg,aes(k_100, SAR))+
  geom_point(aes(color=Bioregions))

cor.test(k_tablePb$k_100, k_tablePb$SAR, method=c("pearson"))

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



# k vs time frame fitting figure

std <- function(x) sd(x)/sqrt(length(x))

fit_figPb<-ggplot(k_tablePb_Sg, aes( Max_Age, k_m2000))+ ggtitle("Decay rate by time frame") + xlab("Time frame (years)") + ylab("Decay rate (yr-1)") +
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
  
  geom_point(data= k_rev,mapping = aes(k_rev$Max_Age, k_rev$k), color='green', shape= 17)+
  
  geom_line(data= fitY1, mapping = aes(timeframe, predict), color='black')+
  
  xlim(0,4000)+
  
  annotate("text", x=1500, y=0.02, label= expression(y == 0.027 * e ** (-0.0018 * 
                                                                              x)))

ggsave(
  plot = fit_figPb,
  path = Folder,
  filename = "fit_plot.jpg",
  units = "cm",
  width = 12,
  height = 7
)



# results distribution

#get coordinates from B dataframe

WM <- map_data("world")

mapa<-k_tablePb[rowSums(is.na(k_tablePb)) != ncol(k_tablePb)-2,]

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






