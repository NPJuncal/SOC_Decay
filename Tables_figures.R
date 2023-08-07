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
DT2 <-merge(DT2, SAR[,c(1,3:6)], by = "ID", all = TRUE)

cor.test(DT2$Av_Mud_25, DT2$Av_13C_25)
plot(DT2$Av_Mud_25, DT2$Av_13C_25)


#### Figure 1 ####

      #Seagrass meadows
      DT2Sg<-subset(DT2, Ecosystem=="Seagrass")
      
      
      cor.test(DT2Sg$Av_C_25, DT2Sg$Av_13C_25, method="spearman")
      plot(DT2Sg$Av_C_25, DT2Sg$Av_13C_25)
      
      #organic carbon  
      
      SC<-ggplot(DT2Sg, aes(C_Gr, Av_C_25)) + ggtitle("Seagrass")+
        geom_boxplot()+
        geom_jitter(color="green4")+
        ylim(-5,50)+
        annotate("text",
                 x = 1:length(table(DT2Sg$C_Gr)),
                 y = -5,
                 label = table(subset(DT2Sg, !is.na(Av_C_25))[,"C_Gr"]),
                 col = "black")+
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      
      pairwise.wilcox.test(DT2Sg$Av_C_25, DT2Sg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)
      #mud  
      SM<-ggplot(DT2Sg, aes(C_Gr, Av_Mud_25)) + 
        geom_boxplot()+
        geom_jitter(color="green4")+
        ylim(-10,110)+
        annotate("text",
                 x = 1:length(table(DT2Sg$C_Gr)),
                 y = -7,
                 label = table(subset(DT2Sg, !is.na(Av_Mud_25))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      
      pairwise.wilcox.test(DT2Sg$Av_Mud_25, DT2Sg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)  
      
      #SAR
      
      SS<-ggplot(DT2Sg, aes(C_Gr, SAR)) + 
        geom_boxplot()+
        geom_jitter(color="green4")+
        ylim(-0.15,1.1)+
        annotate("text",
                 x = 1:length(table(DT2Sg$C_Gr)),
                 y = -0.1,
                 label = table(subset(DT2Sg, !is.na(SAR))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      
      pairwise.wilcox.test(DT2Sg$SAR, DT2Sg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)  
      
      #d13C
      
      d13S<-ggplot(DT2Sg, aes(C_Gr, Av_13C_25)) + 
        geom_boxplot()+
        geom_jitter(color="green4")+
        ylim(-30,-5)+
        annotate("text",
                 x = 1:length(table(DT2Sg$C_Gr)),
                 y = -30,
                 label = table(subset(DT2Sg, !is.na(Av_13C_25))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      
      pairwise.wilcox.test(DT2Sg$Av_13C_25, DT2Sg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)  
      
      
      #Tidal marshes  
      DT2Sm<-subset(DT2, Ecosystem=="Tidal Marsh")
      
      #organic carbon  
      TMC<-ggplot(DT2Sm, aes(C_Gr, Av_C_25)) + ggtitle("Tidal Marsh")+
        geom_boxplot()+
        geom_jitter(color="orange")+
        ylim(-5,50)+
        annotate("text",
                 x = 1:length(table(DT2Sm$C_Gr)),
                 y = -5,
                 label = table(subset(DT2Sm, !is.na(Av_C_25))[,"C_Gr"]),
                 col = "black")+
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      pairwise.wilcox.test(DT2Sm$Av_C_25, DT2Sm$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)
      #mud  
      TMM<-ggplot(DT2Sm, aes(C_Gr, Av_Mud_25)) +
        geom_boxplot()+
        geom_jitter(color="orange")+
        ylim(-10,110)+
        annotate("text",
                 x = 1:length(table(DT2Sm$C_Gr)),
                 y = -7,
                 label = table(subset(DT2Sm, !is.na(Av_Mud_25))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      pairwise.wilcox.test(DT2Sm$Av_Mud_25, DT2Sm$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)  
      
      #SAR 
      TMS<-ggplot(DT2Sm, aes(C_Gr, SAR)) +
        geom_boxplot()+
        geom_jitter(color="orange")+
        ylim(-0.15,1.1)+
        annotate("text",
                 x = 1:length(table(DT2Sm$C_Gr)),
                 y = -0.1,
                 label = table(subset(DT2Sm, !is.na(SAR))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      pairwise.wilcox.test(DT2Sm$SAR, DT2Sm$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05) 
      
      #d13C
      
      d13T<-ggplot(DT2Sm, aes(C_Gr, Av_13C_25)) + 
        geom_boxplot()+
        geom_jitter(color="orange")+
        ylim(-30,-5)+
        annotate("text",
                 x = 1:length(table(DT2Sm$C_Gr)),
                 y = -30,
                 label = table(subset(DT2Sm, !is.na(Av_13C_25))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank())
      
      
      pairwise.wilcox.test(DT2Sm$Av_13C_25, DT2Sm$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)
      
      
      #Mangroves        
      DT2Mg<-subset(DT2, Ecosystem=="Mangrove")
      #organic carbon  
      MGC<-ggplot(DT2Mg, aes(C_Gr, Av_C_25)) + ylab("OC% (Top 25cm)") + ggtitle("Mangrove")+
        geom_boxplot()+
        geom_jitter(color="blue")+
        ylim(-5,50)+
        annotate("text",
                 x = 1:length(table(DT2Mg$C_Gr)),
                 y = -5,
                 label = table(subset(DT2Mg, !is.na(Av_C_25))[,"C_Gr"]),
                 col = "black")+
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank())
      
      pairwise.wilcox.test(DT2Mg$Av_C_25, DT2Mg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05)
      #mud  
      TGM<-ggplot(DT2Mg, aes(C_Gr, Av_Mud_25)) + ylab("Mud% (Top 25cm)") +
        geom_boxplot()+
        geom_jitter(color="blue")+
        ylim(-10,110)+
        annotate("text",
                 x = 1:length(table(DT2Mg$C_Gr)),
                 y = -7,
                 label = table(subset(DT2Mg, !is.na(Av_Mud_25))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank())
      
      pairwise.wilcox.test(DT2Mg$Av_Mud_25, DT2Mg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05) 
      
      #SAR 
      TGS<-ggplot(DT2Mg, aes(C_Gr, SAR)) + ylab("Sediment acc. rate (cm)") +
        geom_boxplot()+
        geom_jitter(color="blue")+
        ylim(-0.15,1.1)+
        annotate("text",
                 x = 1:length(table(DT2Mg$C_Gr)),
                 y = -0.1,
                 label = table(subset(DT2Mg, !is.na(SAR))[,"C_Gr"]),
                 col = "black")+
        theme(axis.title.x = element_blank())
      
      pairwise.wilcox.test(DT2Mg$SAR, DT2Mg$C_Gr,
                           p.adjust.method = "BH") # are significantly different (p < 0.05) 
      
      #d13C
      
      d13M<-ggplot(DT2Mg, aes(C_Gr, Av_13C_25)) + ylab("Average d13C") +
        geom_boxplot()+
        geom_jitter(color="blue")+
        #ylim(-0.15,1.1)+
        #annotate("text",
        #        x = 1:length(table(DT2Mg$C_Gr)),
        #       y = -0.1,
        #      label = table(subset(DT2Mg, !is.na(Av_13C_25))[,"C_Gr"]),
        #     col = "black")+
        theme(axis.title.x = element_blank())
      
      
      #pairwise.wilcox.test(DT2Mg$Av_13C_25, DT2Mg$C_Gr,
      #                    p.adjust.method = "BH") # are significantly different (p < 0.05)
      
      
      todos_CM<-ggpubr::ggarrange(MGC, SC, TMC, TGM, SM, TMM, TGS, SS, TMS,d13M, d13S, d13T,
                                  labels = c("A","B", "C","D" , "E","F", "G","H", "I", "J", "K", "L"),
                                  ncol=3, nrow= 4)
      
      
      ggsave( plot = todos_CM,
              path = Folder,
              filename =  "C_Mud_Gr.jpg",
              units = "cm",
              width = 20,
              height = 25
      )    





#### Table 1 ####

#comparisong among methods and final tables and figures

k_tablePb$G<-"Pb"
k_tableC$G<-"C"

k_tableF<-rbind(k_tablePb,k_tableC)
k_tableF$Ecosystem<-"NA"

for (j in 1:nrow(k_tableF)) {
  
Eco <- substr(k_tableF[j,"ID"], 1, 2)

if (Eco == "Mg") { k_tableF[j, "Ecosystem"] <- "Mangrove"}
if (Eco == "Sg") { k_tableF[j, "Ecosystem"] <- "Seagrass"}
if (Eco == "Sm") { k_tableF[j, "Ecosystem"] <- "Tidal Marsh"}}

k_tableFm<-melt (k_tableF[, c(1:8, 11, 12)], id=c( "ID", "G", "Ecosystem" ))

ggplot(k_tableFm, aes(G, value))+
  geom_boxplot()+
  geom_jitter(aes(color=Ecosystem))+
  ylim(0,0.1)+
  facet_wrap(~variable)

pairwise.wilcox.test(k_tableF$k_2000, k_tableF$G,
                     p.adjust.method = "BH")


# Final table (Pb <500, 14C >500) -----------------------------------------

table_k<-merge(k_tablePb, k_tableC, by = "ID", all = TRUE)
table_k<-table_k[,c(1:5,16:20)]


table_k$Ecosystem<-"NA"

for (j in 1:nrow(table_k)) {
  
  Eco <- substr(table_k[j,"ID"], 1, 2)
  
  if (Eco == "Mg") { table_k[j, "Ecosystem"] <- "Mangrove"}
  if (Eco == "Sg") { table_k[j, "Ecosystem"] <- "Seagrass"}
  if (Eco == "Sm") { table_k[j, "Ecosystem"] <- "Tidal Marsh"}}

colnames(table_k)<-c("ID", "k100", "k150", "k300", "k500", "k1000", "k1500", "k2000", "m2000", "MaxAge", "Ecosystem")
table_k<-table_k[!is.na(table_k$ID),]



#normal distribution and significant differences among ecosystems

shapiro.test(table_k$k150) #normal if pvalue > than 0.05

apply(table_k[,c(2:8)], FUN=shapiro.test, MARGIN = 2)


# differences among species

pairwise.wilcox.test(table_k$k1000, table_k$Ecosystem,
                     p.adjust.method = "BH")
temp<-subset(table_k, !Ecosystem=="Mangrove")
pairwise.wilcox.test(temp$k1000, temp$Ecosystem,
                     p.adjust.method = "BH")

# summary table (manuscript Table 1) -----------------------------------


table_k_Mg<-subset(table_k, Ecosystem=='Mangrove')
table_k_Mg<-sapply(table_k_Mg,FUN=as.numeric)
table_k_Sg<-subset(table_k, Ecosystem=='Seagrass')
table_k_Sg<-sapply(table_k_Sg,FUN=as.numeric)
table_k_Sm<-subset(table_k, Ecosystem=='Tidal Marsh')
table_k_Sm<-sapply(table_k_Sm,FUN=as.numeric)

std <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
cnt <- function(x) sum(!is.na(x))   

sum_table<-as.data.frame(colMeans(table_k_Mg[,c(2:8)], na.rm = TRUE))
sum_table[,2]<-as.data.frame(apply(table_k_Mg[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,3]<-as.data.frame(apply(table_k_Mg[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,4]<-colMeans(table_k_Sg[,c(2:8)], na.rm = TRUE)
sum_table[,5]<-as.data.frame(apply(table_k_Sg[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,6]<-as.data.frame(apply(table_k_Sg[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,7]<-colMeans(table_k_Sm[,c(2:8)], na.rm = TRUE)
sum_table[,8]<-as.data.frame(apply(table_k_Sm[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,9]<-as.data.frame(apply(table_k_Sm[,c(2:8)], FUN=cnt, MARGIN = 2))

sum_table[,10]<-colMeans(table_k[,c(2:8)], na.rm = TRUE)
sum_table[,11]<-as.data.frame(apply(table_k[,c(2:8)], FUN=std, MARGIN = 2))
sum_table[,12]<-as.data.frame(apply(table_k[,c(2:8)], FUN=cnt, MARGIN = 2))

colnames(sum_table)<-c("Mean Mangrove", "SE Mangrove", "n Mangrove","Mean Seagrass", "SE Seagrass", "n Seagrass",
                       "Mean Tidal Marsh", "SE Tidal Marsh", "n Tidal Marsh", "Mean All", "SE All", "n All")

Folder = "Decay2023"
write.csv(sum_table,
          file.path(Folder, "Summar decay.csv"),
          sep = ";",
          dec = ".")


# boxplot by timeframe figure ---------------------------------------------


mtable_k<-melt(table_k[,c(1:8, 11)], id = c("ID","Ecosystem")) 

mtable_k$variable <- as.character(mtable_k$variable)
mtable_k$variable[mtable_k$variable == 'k100'] <- '80-100 yr'
mtable_k$variable[mtable_k$variable == 'k150'] <- '100-150 yr'
mtable_k$variable[mtable_k$variable == 'k300'] <- '150-300 yr'
mtable_k$variable[mtable_k$variable == 'k500'] <- '300-500 yr'
mtable_k$variable[mtable_k$variable == 'k1000'] <- '500-1000 yr'
mtable_k$variable[mtable_k$variable == 'k1500'] <- '1000-1500 yr'
mtable_k$variable[mtable_k$variable == 'k2000'] <- '1500-2000 yr'

mk_table$value<-as.numeric(mk_table$value)

ggplot(transform(mtable_k,
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


mtable_k2<-subset(mtable_k, mtable_k$variable == "100-150 yr" | mtable_k$variable == "500-1000 yr")

box_150_1000<-  ggplot(transform(mtable_k2,
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

f_table_Sg<-as.data.frame(colMeans(table_k[c(7:90),c(2:8)], na.rm = TRUE))
f_table_Sg[,2]<-c(90, 125, 225, 400, 750, 1250, 1750)
f_table_Sg[c(8:13),1]<-na.omit(table_k[,9])
f_table_Sg[c(8:13),2]<-na.omit(table_k[,10])
f_table_Sg[c(14:18), c(1,2)]<-k_revS[,c(2:3)]

#f_table_Sg<-f_table[-9,]
colnames(f_table_Sg)<-c("k","timeframe")

f_table_Tm<-as.data.frame(colMeans(table_k[c(91:115),c(2:6)], na.rm = TRUE))
f_table_Tm[,2]<-c(90, 125, 225, 400, 750)
f_table_Tm[c(6:7), c(1,2)]<-k_revM[,c(2:3)]
colnames(f_table_Tm)<-c("k","timeframe")

f_table_Mg<-as.data.frame(colMeans(table_k[c(1:6),c(2:5)], na.rm = TRUE))
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

fit_fig<-  ggplot(table_k_Sg, aes( MaxAge, m2000))+ ggtitle("Decay rate by time frame") + xlab("Time frame (years)") + ylab("Decay rate (yr-1)") +
  geom_point(color='green4')+
  geom_point(aes(50, mean(na.omit(table_k_Sg$k100))), color='green4')+
  geom_errorbar(aes(50, ymin=mean(na.omit(table_k_Sg$k100))-std(na.omit(table_k_Sg$k100)), ymax=mean(na.omit(table_k_Sg$k100))+std(na.omit(table_k_Sg$k100))), color='green4')+
  geom_point(aes(50, mean(na.omit(table_k_Mg$k100))), color='blue')+
  geom_errorbar(aes(50, ymin=mean(na.omit(table_k_Mg$k100))-std(na.omit(table_k_Mg$k100)), ymax=mean(na.omit(table_k_Mg$k100))+std(na.omit(table_k_Mg$k100))), color='blue')+
  geom_point(aes(50, mean(na.omit(table_k_Sm$k100))), color='orange')+
  geom_errorbar(aes(50, ymin=mean(na.omit(table_k_Sm$k100))-std(na.omit(table_k_Sm$k100)), ymax=mean(na.omit(table_k_Sm$k100))+std(na.omit(table_k_Sm$k100))), color='orange')+
  
  geom_point(aes(125, mean(na.omit(table_k_Sg$k150))), color='green4')+
  geom_errorbar(aes(125, ymin=mean(na.omit(table_k_Sg$k150))-std(na.omit(table_k_Sg$k150)), ymax=mean(na.omit(table_k_Sg$k150))+std(na.omit(table_k_Sg$k150))), color='green4')+
  geom_point(aes(125, mean(na.omit(table_k_Mg$k150))), color='blue')+
  geom_errorbar(aes(125, ymin=mean(na.omit(table_k_Mg$k150))-std(na.omit(table_k_Mg$k150)), ymax=mean(na.omit(table_k_Mg$k150))+std(na.omit(table_k_Mg$k150))), color='blue')+  
  geom_point(aes(125, mean(na.omit(table_k_Sm$k150))), color='orange')+
  geom_errorbar(aes(125, ymin=mean(na.omit(table_k_Sm$k150))-std(na.omit(table_k_Sm$k150)), ymax=mean(na.omit(table_k_Sm$k150))+std(na.omit(table_k_Sm$k150))), color='orange')+  
  
  geom_point(aes(175, mean(na.omit(table_k_Sg$k300))), color='green4')+
  geom_errorbar(aes(175, ymin=mean(na.omit(table_k_Sg$k300))-std(na.omit(table_k_Sg$k300)), ymax=mean(na.omit(table_k_Sg$k300))+std(na.omit(table_k_Sg$k300))), color='green4')+
  geom_point(aes(175, mean(na.omit(table_k_Mg$k300))), color='blue')+
  geom_errorbar(aes(175, ymin=mean(na.omit(table_k_Mg$k300))-std(na.omit(table_k_Mg$k300)), ymax=mean(na.omit(table_k_Mg$k300))+std(na.omit(table_k_Mg$k300))), color='blue')+
  geom_point(aes(175, mean(na.omit(table_k_Sm$k300))), color='orange')+
  geom_errorbar(aes(175, ymin=mean(na.omit(table_k_Sm$k300))-std(na.omit(table_k_Sm$k300)), ymax=mean(na.omit(table_k_Sm$k300))+std(na.omit(table_k_Sm$k300))), color='orange')+
  
  geom_point(aes(400, mean(na.omit(table_k_Sg$k500))), color='green4')+
  geom_errorbar(aes(400, ymin=mean(na.omit(table_k_Sg$k500))-std(na.omit(table_k_Sg$k500)), ymax=mean(na.omit(table_k_Sg$k500))+std(na.omit(table_k_Sg$k500))), color='green4')+
  geom_point(aes(400, mean(na.omit(table_k_Mg$k500))), color='blue')+
  geom_errorbar(aes(400, ymin=mean(na.omit(table_k_Mg$k500))-std(na.omit(table_k_Mg$k500)), ymax=mean(na.omit(table_k_Mg$k500))+std(na.omit(table_k_Mg$k500))), color='blue')+
  geom_point(aes(400, mean(na.omit(table_k_Sm$k500))), color='orange')+
  geom_errorbar(aes(400, ymin=mean(na.omit(table_k_Sm$k500))-std(na.omit(table_k_Sm$k500)), ymax=mean(na.omit(table_k_Sm$k500))+std(na.omit(table_k_Sm$k500))), color='orange')+
  
  geom_point(aes(750, mean(na.omit(table_k_Sg$k1000))), color='green4')+
  geom_errorbar(aes(750, ymin=mean(na.omit(table_k_Sg$k1000))-std(na.omit(table_k_Sg$k1000)), ymax=mean(na.omit(table_k_Sg$k1000))+std(na.omit(table_k_Sg$k1000))), color='green4')+
  geom_point(aes(750, mean(na.omit(table_k_Mg$k1000))), color='blue')+
  geom_errorbar(aes(750, ymin=mean(na.omit(table_k_Mg$k1000))-std(na.omit(table_k_Mg$k1000)), ymax=mean(na.omit(table_k_Mg$k1000))+std(na.omit(table_k_Mg$k1000))), color='blue')+
  geom_point(aes(750, mean(na.omit(table_k_Sm$k1000))), color='orange')+
  geom_errorbar(aes(750, ymin=mean(na.omit(table_k_Sm$k1000))-std(na.omit(table_k_Sm$k1000)), ymax=mean(na.omit(table_k_Sm$k1000))+std(na.omit(table_k_Sm$k1000))), color='orange')+
  
  geom_point(aes(1250, mean(na.omit(table_k_Sg$k1500))), color='green4')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(table_k_Sg$k1500))-std(na.omit(table_k_Sg$k1500)), ymax=mean(na.omit(table_k_Sg$k1500))+std(na.omit(table_k_Sg$k1500))), color='green4')+
  geom_point(aes(1250, mean(na.omit(table_k_Mg$k1500))), color='blue')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(table_k_Mg$k1500))-std(na.omit(table_k_Mg$k1500)), ymax=mean(na.omit(table_k_Mg$k1500))+std(na.omit(table_k_Mg$k1500))), color='blue')+
  geom_point(aes(1250, mean(na.omit(table_k_Sm$k1500))), color='orange')+
  geom_errorbar(aes(1250, ymin=mean(na.omit(table_k_Sm$k1500))-std(na.omit(table_k_Sm$k1500)), ymax=mean(na.omit(table_k_Sm$k1500))+std(na.omit(table_k_Sm$k1500))), color='orange')+
  
  geom_point(aes(1750, mean(na.omit(table_k_Sg$k2000))), color='green4')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(table_k_Sg$k2000))-std(na.omit(table_k_Sg$k2000)), ymax=mean(na.omit(table_k_Sg$k2000))+std(na.omit(table_k_Sg$k2000))), color='green4')+
  geom_point(aes(1750, mean(na.omit(table_k_Mg$k2000))), color='blue')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(table_k_Mg$k2000))-std(na.omit(table_k_Mg$k2000)), ymax=mean(na.omit(table_k_Mg$k2000))+std(na.omit(table_k_Mg$k2000))), color='blue')+
  geom_point(aes(1750, mean(na.omit(table_k_Sm$k2000))), color='orange')+
  geom_errorbar(aes(1750, ymin=mean(na.omit(table_k_Sm$k2000))-std(na.omit(table_k_Sm$k2000)), ymax=mean(na.omit(table_k_Sm$k2000))+std(na.omit(table_k_Sm$k2000))), color='orange')+
  
  geom_point(data= k_revS,mapping = aes(k_revS$Max_Age, k_revS$k), color='green', shape= 17)+
  geom_point(data= k_revM,mapping = aes(k_revM$Max_Age, k_revM$k), color='orange4', shape= 17)+
  
  geom_line(data= fitSg, mapping = aes(timeframe, predict), color='green4')+
  geom_line(data= fitTm, mapping = aes(timeframe, predict), color='orange')+
  geom_line(data= fitMg, mapping = aes(timeframe, predict), color='blue')+
  
  xlim(0,4000)+
  
  annotate("text", x=1500, y=0.025, color= "blue",  size = 5, label= expression(y == 0.028 * e ** (-0.004 * 
                                                                                                     x)))+
  annotate("text", x=1500, y=0.02, color= "green4",size = 5,label= expression(y == 0.037 * e ** (-0.003 * 
                                                                                                   x)))+
  annotate("text", x=1500, y=0.015, color= "orange",size = 5,label= expression(y == 0.023 * e ** (-0.003 * 
                                                                                                    x)))


ggsave(
  plot = fit_fig,
  path = Folder,
  filename = "fit_plot.jpg",
  units = "cm",
  width = 12,
  height = 7
)



# correlation decay rate, OC, Mud and d13C ------------------------------------

# first we estimate the average OC, mud, d13C and SAR content for the studied time frame
# we use the max age of each core

x<-split(TPb, TPb$Core)

Sum150_a<-fit_150Pb
Sum150_a$OC<-"NA"
Sum150_a$Mud<-"NA"
Sum150_a$d13C<-"NA"
Sum150_a$SAR<-"NA"

for (i in 1:nrow(Sum150_a)) {
  
  data<-as.data.frame(x[Sum150_a[i,"ID"]])
  colnames(data)<-colnames(TPb)
  data_a<-subset(data, data$FAge<=Sum150_a[i, "Max.Age"])
  Sum150_a[i, "OC"]<-mean(data_a$Corg, na.rm=TRUE)
  Sum150_a[i, "Mud"]<-mean(data_a$Mud, na.rm=TRUE)
  Sum150_a[i, "d13C"]<-mean(data_a$d13C, na.rm=TRUE)
  Sum150_a[i, "SAR"]<-max(data_a$Max.Depth)/max(data_a$FAge)}

Sum150_a$Mud<-as.numeric(Sum150_a$Mud)
Sum150_a$OC<-as.numeric(Sum150_a$OC)
Sum150_a$d13C<-as.numeric(Sum150_a$d13C)
Sum150_a$SAR<-as.numeric(Sum150_a$SAR)

Sum1000_a<-fit_1000C
Sum1000_a$OC<-"NA"
Sum1000_a$Mud<-"NA"
Sum1000_a$d13C<-"NA"
Sum1000_a$SAR<-"NA"

x<-split(TC, TC$Core)

for (i in 1:nrow(Sum1000_a)) {
  
  data<-as.data.frame(x[Sum1000_a[i,"ID"]])
  colnames(data)<-colnames(TC)
  data_a<-subset(data, data$FAge<=Sum1000_a[i, "Max.Age"])
  Sum1000_a[i, "OC"]<-mean(data_a$Corg, na.rm=TRUE)
  Sum1000_a[i, "Mud"]<-mean(data_a$Mud, na.rm=TRUE)
  Sum1000_a[i, "d13C"]<-mean(data_a$d13C, na.rm=TRUE)
  Sum1000_a[i, "SAR"]<-max(data_a$Max.Depth)/max(data_a$FAge)}

Sum1000_a$Mud<-as.numeric(Sum1000_a$Mud)
Sum1000_a$OC<-as.numeric(Sum1000_a$OC)
Sum1000_a$d13C<-as.numeric(Sum1000_a$d13C)
Sum1000_a$SAR<-as.numeric(Sum1000_a$SAR)


shapiro.test(Sum150_a$k) #normal if pvalue > than 0.05
shapiro.test(as.numeric(Sum150_a$OC))
shapiro.test(as.numeric(Sum150_a$Mud))
shapiro.test(as.numeric(Sum150_a$d13C))
shapiro.test(as.numeric(Sum150_a$SAR))
plot(Sum150_a$k, Sum150_a$Mud)

cor.test(as.numeric(Sum150_a$k), as.numeric(Sum150_a$Mud), method=c("spearman"))
cor.test(as.numeric(Sum150_a$k), as.numeric(Sum150_a$OC), method=c("spearman"))
cor.test(as.numeric(Sum150_a$k), as.numeric(Sum150_a$d13C), method=c("spearman"))
cor.test(as.numeric(Sum150_a$k), as.numeric(Sum150_a$SAR), method=c("spearman"))



ggplot(Sum150_a, aes(k, OC))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c('blue', 'green4', "orange"))

ggplot(Sum150_a, aes(k, d13C))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c('blue', 'green4', "orange"))

ggplot(Sum150_a, aes(k, SAR))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c('blue', 'green4', "orange"))



shapiro.test(Sum1000_a$k) #normal if pvalue > than 0.05
shapiro.test(as.numeric(Sum1000_a$Mud))
shapiro.test(as.numeric(Sum1000_a$d13C))
shapiro.test(as.numeric(Sum1000_a$SAR))
cor.test(as.numeric(Sum1000_a$k), as.numeric(Sum1000_a$Mud), method=c("spearman"))
cor.test(as.numeric(Sum1000_a$k), as.numeric(Sum1000_a$OC), method=c("spearman"))
cor.test(as.numeric(Sum1000_a$k), as.numeric(Sum1000_a$d13C), method=c("spearman"))
cor.test(as.numeric(Sum1000_a$k), as.numeric(Sum1000_a$SAR), method=c("spearman"))

ltm<-ggplot(Sum1000_a, aes(k, Mud))+ ylab("Mud % (<0.063 mm)") +
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c( 'green4', "orange"))

ggplot(Sum1000_a, aes(k, OC))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c( 'green4', "orange"))

ggplot(Sum1000_a, aes(k, d13C))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c( 'green4', "orange"))

ggplot(Sum1000_a, aes(k, SAR))+
  geom_point(aes(color=Ecosystem))+
  scale_color_manual(values=c('green4', "orange"))

# figure 3
stm<-ggplot(Sum150_a, aes(k, Mud))+ ylab("Mud % (<0.063 mm)") + xlab("100-150 yr") +
  geom_point(aes(color=Ecosystem))+
  xlim(0, 0.04)+ ylim (0,100)+
  scale_color_manual(values=c('blue', 'green4', "orange"))

ltm<-ggplot(Sum1000_a, aes(k, Mud))+ ylab("Mud % (<0.063 mm)") + xlab("500-1000 yr") +
  geom_point(aes(color=Ecosystem))+
  xlim(0, 0.04)+ ylim (0,100)+
  scale_color_manual(values=c( 'green4', "orange"))


mud_decay<-grid.arrange(stm, ltm, top = "Dacay rates (yr-1)")


ggsave(
  plot = mud_decay,
  path = Folder,
  filename = "mud_decay.jpg",
  units = "cm",
  width = 9,
  height = 10
)



# boxplots and wilcox -------------------------------------------------












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











# # results distribution (Fig 4) --------------------------------------------------


#get coordinates from B dataframe

WM <- map_data("world")

mapa<-k_tablePb[rowSums(is.na(k_tablePb)) != ncol(k_tablePb)-2,]

CoordR<-subset(B, B$Core %in% mapa$ID==TRUE)
CoordR$Lat<-as.numeric(CoordR$Lat)
CoordR$Long<-as.numeric(CoordR$Long)

    global2<-CoordR %>%
      ggplot() + ggtitle("Estimated k distribution") + xlab("Longitude") + ylab("Latitude") +
      geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
      #geom_point(aes(x = Long, y = Lat))+
      geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
                 pch = 21,
                 size = 1.8) +
      coord_sf(xlim = c(-140, 150), ylim = c(-40, 75)) +
      scale_fill_manual(values = c( "blue", "green","orange")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    nam2<-CoordR %>%
      ggplot() + xlab("Longitude") + ylab("Latitude") +
      geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
      #geom_point(aes(x = Long, y = Lat))+
      geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
                 pch = 21,
                 size = 1.8) +
      coord_sf(xlim = c(-150, -50), ylim = c(0, 70)) +
      scale_fill_manual(values = c( "blue", "green","orange")) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none")
    
    eu2<-CoordR %>%
      ggplot() + xlab("Longitude") + ylab("Latitude") +
      geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
      #geom_point(aes(x = Long, y = Lat))+
      geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
                 pch = 21,
                 size = 1.8) +
      coord_sf(xlim = c(-10, 50), ylim = c(20, 70)) +
      scale_fill_manual(values = c( "blue", "green","orange")) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none")
    
    aus2<-CoordR %>%
      ggplot() + xlab("Longitude") + ylab("Latitude") +
      geom_polygon(data = WM, aes(x = long, y = lat, group = group)) +
      #geom_point(aes(x = Long, y = Lat))+
      geom_point(aes(x = Long, y = Lat,  fill = Ecosystem),
                 pch = 21,
                 size = 1.8) +
      coord_sf(xlim = c(110, 152), ylim = c(-45, -10))  +
      scale_fill_manual(values = c( "blue", "green","orange")) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none")



    ssp2<-grid.arrange(global2, nam2, eu2, aus2, 
                      layout_matrix = rbind(c(1, 1, 1),
                                            c(2, 3, 4)))
    
    
    ggsave(
      plot = ssp2,
      path = Folder,
      filename =  "estimated k all map.jpg",
      units = "cm",
      width = 20,
      height = 15
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



# from time to decay rate -------------------------------------------------

# degradation of the first meter
#seagrass
0.037 * exp(-0.003 * 1953.206)

0.037 * exp(-0.003 * 1504.5)
0.037 * exp(-0.003 * 2401.9)

#mangrove
0.028 * exp(-0.004 * 1953.206)

0.028 * exp(-0.004 * 1504.5)
0.028 * exp(-0.004 * 2401.9)

#tidal marshes
0.023 * exp(-0.003 * 1953.206)

0.023 * exp(-0.003 * 1504.5)
0.023 * exp(-0.003 * 2401.9)
