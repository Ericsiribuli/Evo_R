# day3信号分析-----

TGY_day3 <- read.table('freshwt_TGY_day3_fitness_single_8-2.txt',sep = ' ')
names(TGY_day3) <- c('mut','fitness_TGY')
TGY_day3$fitness_TGY <- TGY_day3$fitness_TGY ** (1/43)

TGS_day3 <- read.table('freshwt_TGS_day3_fitness_single_8-2.txt',sep = ' ')
names(TGS_day3) <- c('mut','fitness_TGS')
TGS_day3$fitness_TGS <- TGS_day3$fitness_TGS ** (1/38)

AGY_day3 <- read.table('freshwt_AGY_day3_fitness_single_8-2.txt',sep = ' ')
names(AGY_day3) <- c('mut','fitness_AGY')
AGY_day3$fitness_AGY <- AGY_day3$fitness_AGY ** (1/43)

AGS_day3 <- read.table('freshwt_AGS_day3_fitness_single_8-2.txt',sep = ' ')
names(AGS_day3) <- c('mut','fitness_AGS')
AGS_day3$fitness_AGS <- AGS_day3$fitness_AGS ** (1/38)


TGY_TGS_day3 <- merge(TGY_day3,TGS_day3,by='mut',all = F)
wilcox.test(TGY_TGS_day3$fitness_TGY,TGY_TGS_day3$fitness_TGS,paired = T)
median(TGY_TGS_day3$fitness_TGY)
median(TGY_TGS_day3$fitness_TGS)

TGY_AGY_day3 <- merge(TGY_day3,AGY_day3,by='mut',all = F)
wilcox.test(TGY_AGY_day3$fitness_TGY,TGY_AGY_day3$fitness_AGY,paired = T)
median(TGY_AGY_day3$fitness_TGY)
median(TGY_AGY_day3$fitness_AGY)

AGY_AGS_day3 <- merge(AGY_day3,AGS_day3,by='mut',all = F)
wilcox.test(AGY_AGS_day3$fitness_AGY,AGY_AGS_day3$fitness_AGS,paired = T)
median(AGY_AGS_day3$fitness_AGY)
median(AGY_AGS_day3$fitness_AGS)

TGS_AGS_day3 <- merge(TGS_day3,AGS_day3,by='mut',all = F)
wilcox.test(TGS_AGS_day3$fitness_TGS,TGS_AGS_day3$fitness_AGS,paired = T)
median(TGS_AGS_day3$fitness_TGS)
median(TGS_AGS_day3$fitness_AGS)

ggplot(TGY_TGS_AGY_AGS_day3,aes(x =TGY_TGS_AGY_AGS_day3$fitness_TGY,TGY_TGS_AGY_AGS_day3$fitness_AGY)) +
  geom_point(size=1) + xlim(0.95,1.05) + ylim(0.95,1.05) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TGY") + ylab("s_AGY") +
  style.print()

ggplot(TGY_TGS_AGY_AGS_day3, aes(fitness_TGY, fitness_TGS)) +geom_boxplot() +
  ylim(0.95,1.05) +
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=13), axis.title.y=element_text(size=13),legend.text=element_text(size=15)) +
  labs(title = "TGY AGY samples effect size",subtitle = expression(paste(" P =  ",1.6 %*% 10^-3)))+
  theme(plot.title = element_text(size = 20)) 


TUY_day3 <- read.table('freshwt_TUY_day3_fitness_single_8-2.txt',sep = ' ')
names(TUY_day3) <- c('mut','fitness_TUY')
TUY_day3$fitness_TUY <- TUY_day3$fitness_TUY ** (1/43)
TUY_day3$group <-'TUY'

AUY_day3 <- read.table('freshwt_AUY_day3_fitness_single_8-2.txt',sep = ' ')
names(AUY_day3) <- c('mut','fitness_AUY')
AUY_day3$fitness_AUY <- AUY_day3$fitness_AUY ** (1/43)
AUY_day3$group <-'AUY'

TUS_day3 <- read.table('freshwt_TUS_day3_fitness_single_8-2.txt',sep = ' ')
names(TUS_day3) <- c('mut','fitness_TUS')
TUS_day3$fitness_TUS <- TUS_day3$fitness_TUS ** (1/38)
TUS_day3$group <-'TUS'

AUS_day3 <- read.table('freshwt_AUS_day3_fitness_single_8-2.txt',sep = ' ')
names(AUS_day3) <- c('mut','fitness_AUS')
AUS_day3$fitness_AUS <- AUS_day3$fitness_AUS ** (1/38)
AUS_day3$group <-'AUS'

TUU_day3 <- read.table('freshwt_TUU_day3_fitness_single_8-2.txt',sep = ' ')
names(TUU_day3) <- c('mut','fitness_TUU')
TUU_day3$fitness_TUU <- TUU_day3$fitness_TUU ** (1/38)
TUU_day3$group <-'TUU'

AUU_day3 <- read.table('freshwt_AUU_day3_fitness_single_8-2.txt',sep = ' ')
names(AUU_day3) <- c('mut','fitness_AUU')
AUU_day3$fitness_AUU <- AUU_day3$fitness_AUU ** (1/38)
AUU_day3$group <-'AUU'

TUY_AUY_day3 <- merge(TUY_day3,AUY_day3,by='mut',all = F)
wilcox.test(TUY_AUY_day3$fitness_TUY,TUY_AUY_day3$fitness_AUY,paired = T)
median(TUY_AUY_day3$fitness_TUY)
median(TUY_AUY_day3$fitness_AUY)

TUS_AUS_day3 <- merge(TUS_day3,AUS_day3,by='mut',all = F)
wilcox.test(TUS_AUS_day3$fitness_TUS,TUS_AUS_day3$fitness_AUS,paired = T)
median(TUS_AUS_day3$fitness_TUS)
median(TUS_AUS_day3$fitness_AUS)

TUS_TUU_day3 <- merge(TUS_day3,TUU_day3,by='mut',all = F)
wilcox.test(TUS_TUU_day3$fitness_TUS,TUS_TUU_day3$fitness_TUU,paired = T)
median(TUS_TUU_day3$fitness_TUS)
median(TUS_TUU_day3$fitness_TUU)

AUS_AUU_day3 <- merge(AUS_day3,AUU_day3,by='mut',all = F)
wilcox.test(AUS_AUU_day3$fitness_AUS,AUS_AUU_day3$fitness_AUU,paired = T)
median(AUS_AUU_day3$fitness_AUS)
median(AUS_AUU_day3$fitness_AUU)

TGY_TUY_day3 <- merge(TGY_day3,TUY_day3,by='mut',all = F)
wilcox.test(TGY_TUY_day3$fitness_TGY,TGY_TUY_day3$fitness_TUY,paired = T)
median(TGY_TUY_day3$fitness_TGY)
median(TGY_TUY_day3$fitness_TUY)

TGS_TUS_day3 <- merge(TGS_day3,TUS_day3,by='mut',all = F)
wilcox.test(TGS_TUS_day3$fitness_TGS,TGS_TUS_day3$fitness_TUS,paired = T)
median(TGS_TUS_day3$fitness_TGS)
median(TGS_TUS_day3$fitness_TUS)


TUS_AUS_day3 <- rbind(TUS_day3,AUS_day3)
ggplot(TUS_AUS_day3, aes(group, fitness)) +geom_boxplot() +
  ylim(0.98,1.02) +
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=13), axis.title.y=element_text(size=13),legend.text=element_text(size=15)) +
  labs(title = "TUS AUS samples effect size",subtitle = expression(paste(" P =  ",1.6 %*% 10^-3)))+
  theme(plot.title = element_text(size = 20)) 

TUU_TUS_day3 <- rbind(TUU_day3,TUS_day3)
ggplot(TUU_TUS_day3, aes(group, fitness)) +geom_boxplot() +
  ylim(0.98,1.02) +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) +
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=13), axis.title.y=element_text(size=13),legend.text=element_text(size=15)) +
  labs(title = "TUU TUS samples effect size",subtitle = expression(paste(" P =  ",1.6 %*% 10^-3)))+
  theme(plot.title = element_text(size = 20)) 

TUU_TUS_day3_mer <- merge(TUU_day3,TUS_day3,by='mut',all = T)
ggplot(TUU_TUS_day3_mer,aes(x =TUU_TUS_day3_mer$fitness_TUS,TUU_TUS_day3_mer$fitness_TUU)) +
  geom_point(size=1) + xlim(0.93,1.07) + ylim(0.93,1.07) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TUS") + ylab("s_TUU") +
  style.print()

AUU_AUS_day3 <- rbind(AUU_day3,AUS_day3)
ggplot(AUU_AUS_day3, aes(group, fitness)) +geom_boxplot() +
  ylim(0.98,1.02) +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) +
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=13), axis.title.y=element_text(size=13),legend.text=element_text(size=15)) +
  labs(title = "TUU TUS samples effect size",subtitle = expression(paste(" P =  ",1.6 %*% 10^-3)))+
  theme(plot.title = element_text(size = 20))


# GFP的RSA、疏水亲水参数导入-----

setwd("~/misinteration")
hyb <- read.csv("AA_hyb.csv",header = F, stringsAsFactors = F)
hyb <- hyb[,-1]

gfp_ACC <- read.csv('gfp_pos_ACC.csv',header = T, stringsAsFactors = F)
names(gfp_ACC) <-c('AA_pos','AA_bef','ACC')
Surface_area <- read.csv("surface_area_AA.csv",header = T)
names(Surface_area) <- c("AA_bef","sur_area")
gfp_RSA<- merge(gfp_ACC,Surface_area,by='AA_bef',all=T)
w<-which(gfp_RSA$AA_pos=="65")
gfp_RSA <- gfp_RSA[-w,]
gfp_RSA$ACC <-as.numeric(gfp_RSA$ACC)
gfp_RSA$RSA <- gfp_RSA$ACC/gfp_RSA$sur_area

names(hyb) <- c("AA_bef","hyb_bef")
gfp_AA_bef_hyb <- merge(gfp_RSA,hyb,by='AA_bef',all = T)

gfp_AA_mutafter <- read.csv('gfp_AA_after.csv',header = F)
names(gfp_AA_mutafter) <-c('mut','AA_pos','AA_aft')

names(hyb) <- c("AA_aft","hyb_aft")
gfp_AA_mutafter_hyb<- merge(gfp_AA_mutafter,hyb,by='AA_aft',all = T)

gfp_AAba_RSA <- merge(gfp_AA_bef_hyb,gfp_AA_mutafter_hyb,by='AA_pos',all = T)
gfp_AAba_RSA <- na.omit(gfp_AAba_RSA)
gfp_AAba_RSA$hyb_t <- gfp_AAba_RSA$hyb_aft - gfp_AAba_RSA$hyb_bef


# TGY misinteraction misfolding IDR-----

TGY_day3$mut <- as.character(TGY_day3$mut)

TGY_day3$AA_pos <- 
  TGY_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

TGY_day3_RSA <- merge(TGY_day3,gfp_AAba_RSA,by=c('AA_pos','mut'),all = T)
TGY_day3_RSA <- na.omit(TGY_day3_RSA)


TGY_day3_geno <- TGY_day3_RSA[,c(1,2,3,7,11)]

#TGY_test1 <- TGY_geno[,c(2,5,11,12,13)]
#TGY_day3_geno <- merge(TGY_day3,TGY_test1,by='mut',all = T)
#TGY_day3_geno <- na.omit(TGY_day3_geno)

TGY_day3_surfaced <- filter(TGY_day3_geno, RSA > 0.2)
TGY_day3_surfaced <- TGY_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TGY=mean(fitness_TGY))

ggplot(TGY_day3_surfaced,aes(TGY_day3_surfaced$hyb_t,TGY_day3_surfaced$fitness_TGY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TGY diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(TGY_day3_surfaced$hyb_t,TGY_day3_surfaced$fitness_TGY,method = "s")

TGY_day3_surfaced$hyb_RSA <- TGY_day3_surfaced$hyb_t * TGY_day3_surfaced$RSA
ggplot(TGY_day3_surfaced,aes(TGY_day3_surfaced$hyb_RSA,TGY_day3_surfaced$fitness_TGY)) +
  geom_point() + ylim(0.99,1.02) + xlim(-4,6.5) +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TGY_day3_surfaced$hyb_RSA,TGY_day3_surfaced$fitness_TGY,method = "s")


# misfolding 

mut_misfolding <- read.csv('mut_misfolding.csv',header = F)
names(mut_misfolding) <- c('mut','misfolding')

TGY_day3_misfolding <- merge(TGY_day3_geno,mut_misfolding,by='mut',all = T)
TGY_day3_misfolding <- na.omit(TGY_day3_misfolding)
TGY_day3_buried_misfolding <- filter(TGY_day3_misfolding, RSA<0.2)

ggplot(TGY_day3_buried_misfolding,aes(TGY_day3_buried_misfolding$misfolding,TGY_day3_buried_misfolding$fitness_TGY)) +
  geom_point() +
  annotate('text', label = 'rho=-0.112  P=0.0003', x = 350, y = 1.01, size=7, colour ='#000000')+
  geom_smooth(method = 'lm') + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(TGY_day3_buried_misfolding$misfolding,TGY_day3_buried_misfolding$fitness_TGY,method = "s")


# IDR 

mut_disorder <- read.csv('mut_disorder.csv',header = F)
names(mut_disorder) <- c('mut','disorder_length')

TGY_day3_sur_diso <- merge(TGY_day3_geno,mut_disorder,by='mut',all = T)
TGY_day3_sur_diso <- na.omit(TGY_day3_sur_diso)

ggplot(TGY_day3_sur_diso,aes(TGY_day3_sur_diso$disorder_length,TGY_day3_sur_diso$fitness_TGY)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_day3_sur_diso$disorder_length,TGY_day3_sur_diso$fitness_TGY,method = "s")

TGY_day3_sur_mean_disor <- TGY_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_TGY = mean(fitness_TGY))

ggplot(TGY_day3_sur_mean_disor,aes(TGY_day3_sur_mean_disor$disorder_length,TGY_day3_sur_mean_disor$fitness_TGY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TGY IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TGY_day3_sur_mean_disor$disorder_length,TGY_day3_sur_mean_disor$fitness_TGY,method = "s")



# AGY misinteraction misfolding IDR-----
# misinteraction

AGY_day3$mut <- as.character(AGY_day3$mut)

AGY_day3$AA_pos <- 
  AGY_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

AGY_day3_RSA <- merge(AGY_day3,gfp_AAba_RSA,by=c('AA_pos','mut'),all = T)
AGY_day3_RSA <- na.omit(AGY_day3_RSA)

AGY_day3_geno <- AGY_day3_RSA[,c(1,2,3,7,11)]

AGY_day3_surfaced <- filter(AGY_day3_geno, RSA > 0.2)
AGY_day3_surfaced <- AGY_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AGY=mean(fitness_AGY))

ggplot(AGY_day3_surfaced,aes(AGY_day3_surfaced$hyb_t,AGY_day3_surfaced$fitness_AGY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "AGY diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(AGY_day3_surfaced$hyb_t,AGY_day3_surfaced$fitness_AGY,method = "s")

AGY_day3_surfaced$hyb_RSA <- AGY_day3_surfaced$hyb_t * AGY_day3_surfaced$RSA
ggplot(AGY_day3_surfaced,aes(AGY_day3_surfaced$hyb_RSA,AGY_day3_surfaced$fitness_AGY)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(AGY_day3_surfaced$hyb_RSA,AGY_day3_surfaced$fitness_AGY,method = "s")

# misfolding 

mut_misfolding <- read.csv('mut_misfolding.csv',header = F)
names(mut_misfolding) <- c('mut','misfolding')

AGY_day3_misfolding <- merge(AGY_day3_geno,mut_misfolding,by='mut',all = T)
AGY_day3_misfolding <- na.omit(AGY_day3_misfolding)
AGY_day3_buried_misfolding <- filter(AGY_day3_misfolding, RSA<0.2)

ggplot(AGY_day3_buried_misfolding,aes(AGY_day3_buried_misfolding$misfolding,AGY_day3_buried_misfolding$fitness_AGY)) +
  geom_point() +
  annotate('text', label = 'rho=0.069  P= 0.026', x = 350, y = 1.01, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(AGY_day3_buried_misfolding$misfolding,AGY_day3_buried_misfolding$fitness_AGY,method = "s")

# IDR 

mut_disorder <- read.csv('mut_disorder.csv',header = F)
names(mut_disorder) <- c('mut','disorder_length')

AGY_day3_sur_diso <- merge(AGY_day3_geno,mut_disorder,by='mut',all = T)
AGY_day3_sur_diso <- na.omit(AGY_day3_sur_diso)

ggplot(AGY_day3_sur_diso,aes(AGY_day3_sur_diso$disorder_length,AGY_day3_sur_diso$fitness_AGY)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(AGY_day3_sur_diso$disorder_length,AGY_day3_sur_diso$fitness_AGY,method = "s")

AGY_day3_sur_mean_disor <- AGY_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_AGY = mean(fitness_AGY))

ggplot(AGY_day3_sur_mean_disor,aes(AGY_day3_sur_mean_disor$disorder_length,AGY_day3_sur_mean_disor$fitness_AGY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "AGY IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(AGY_day3_sur_mean_disor$disorder_length,AGY_day3_sur_mean_disor$fitness_AGY,method = "s")



# TGS misinteraction misfolding IDR-----
# misinteraction
TGS_day3$mut <- as.character(TGS_day3$mut)

TGS_day3$AA_pos <- 
  TGS_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

TGS_day3_RSA <- merge(TGS_day3,gfp_AAba_RSA,by=c('AA_pos','mut'),all = T)
TGS_day3_RSA <- na.omit(TGS_day3_RSA)

TGS_day3_geno <- TGS_day3_RSA[,c(1,2,3,7,11)]

TGS_day3_surfaced <- filter(TGS_day3_geno, RSA > 0.2)
TGS_day3_surfaced <- TGS_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TGS=mean(fitness_TGS))

ggplot(TGS_day3_surfaced,aes(TGS_day3_surfaced$hyb_t,TGS_day3_surfaced$fitness_TGS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TGS diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(TGS_day3_surfaced$hyb_t,TGS_day3_surfaced$fitness_TGS,method = "s")

TGS_day3_surfaced$hyb_RSA <- TGS_day3_surfaced$hyb_t * TGS_day3_surfaced$RSA
ggplot(TGS_day3_surfaced,aes(TGS_day3_surfaced$hyb_RSA,TGS_day3_surfaced$fitness_TGS)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TGS_day3_surfaced$hyb_RSA,TGS_day3_surfaced$fitness_TGS,method = "s")


# misfolding 

mut_misfolding <- read.csv('mut_misfolding.csv',header = F)
names(mut_misfolding) <- c('mut','misfolding')

TGS_day3_misfolding <- merge(TGS_day3_geno,mut_misfolding,by='mut',all = T)
TGS_day3_misfolding <- na.omit(TGS_day3_misfolding)
TGS_day3_buried_misfolding <- filter(TGS_day3_misfolding, RSA<0.2)

ggplot(TGS_day3_buried_misfolding,aes(TGS_day3_buried_misfolding$misfolding,TGS_day3_buried_misfolding$fitness_TGS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(TGS_day3_buried_misfolding$misfolding,TGS_day3_buried_misfolding$fitness_TGS,method = "s")

# IDR 

mut_disorder <- read.csv('mut_disorder.csv',header = F)
names(mut_disorder) <- c('mut','disorder_length')

TGS_day3_sur_diso <- merge(TGS_day3_geno,mut_disorder,by='mut',all = T)
TGS_day3_sur_diso <- na.omit(TGS_day3_sur_diso)

ggplot(TGS_day3_sur_diso,aes(TGS_day3_sur_diso$disorder_length,TGS_day3_sur_diso$fitness_TGS)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGS_day3_sur_diso$disorder_length,TGS_day3_sur_diso$fitness_TGS,method = "s")

TGS_day3_sur_mean_disor <- TGS_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_TGS = mean(fitness_TGS))

ggplot(TGS_day3_sur_mean_disor,aes(TGS_day3_sur_mean_disor$disorder_length,TGS_day3_sur_mean_disor$fitness_TGS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TGS IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TGS_day3_sur_mean_disor$disorder_length,TGS_day3_sur_mean_disor$fitness_TGS,method = "s")


# AGS misinteraction misfolding IDR-----
# misinteraction
AGS_day3$mut <- as.character(AGS_day3$mut)

AGS_day3$AA_pos <- 
  AGS_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

AGS_day3_RSA <- merge(AGS_day3,gfp_AAba_RSA,by=c('AA_pos','mut'),all = T)
AGS_day3_RSA <- na.omit(AGS_day3_RSA)

AGS_day3_geno <- AGS_day3_RSA[,c(1,2,3,7,11)]

AGS_day3_surfaced <- filter(AGS_day3_geno, RSA > 0.2)
AGS_day3_surfaced <- AGS_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AGS=mean(fitness_AGS))


ggplot(AGS_day3_surfaced,aes(AGS_day3_surfaced$hyb_t,AGS_day3_surfaced$fitness_AGS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "AGS diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(AGS_day3_surfaced$hyb_t,AGS_day3_surfaced$fitness_AGS,method = "s")

AGS_day3_surfaced$hyb_RSA <- AGS_day3_surfaced$hyb_t * AGS_day3_surfaced$RSA
ggplot(AGS_day3_surfaced,aes(AGS_day3_surfaced$hyb_RSA,AGS_day3_surfaced$fitness_AGS)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(AGS_day3_surfaced$hyb_RSA,AGS_day3_surfaced$fitness_AGS,method = "s")


# misfolding 

mut_misfolding <- read.csv('mut_misfolding.csv',header = F)
names(mut_misfolding) <- c('mut','misfolding')

AGS_day3_misfolding <- merge(AGS_day3_geno,mut_misfolding,by='mut',all = T)
AGS_day3_misfolding <- na.omit(AGS_day3_misfolding)
AGS_day3_buried_misfolding <- filter(AGS_day3_misfolding, RSA<0.2)

ggplot(AGS_day3_buried_misfolding,aes(AGS_day3_buried_misfolding$misfolding,AGS_day3_buried_misfolding$fitness_AGS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(AGS_day3_buried_misfolding$misfolding,AGS_day3_buried_misfolding$fitness_AGS,method = "s")

# IDR 

mut_disorder <- read.csv('mut_disorder.csv',header = F)
names(mut_disorder) <- c('mut','disorder_length')

AGS_day3_sur_diso <- merge(AGS_day3_geno,mut_disorder,by='mut',all = T)
AGS_day3_sur_diso <- na.omit(AGS_day3_sur_diso)

ggplot(AGS_day3_sur_diso,aes(AGS_day3_sur_diso$disorder_length,AGS_day3_sur_diso$fitness_AGS)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(AGS_day3_sur_diso$disorder_length,AGS_day3_sur_diso$fitness_AGS,method = "s")

AGS_day3_sur_mean_disor <- AGS_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_AGS = mean(fitness_AGS))

ggplot(AGS_day3_sur_mean_disor,aes(AGS_day3_sur_mean_disor$disorder_length,AGS_day3_sur_mean_disor$fitness_AGS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "AGS IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(AGS_day3_sur_mean_disor$disorder_length,AGS_day3_sur_mean_disor$fitness_AGS,method = "s")


# URA3的RSA、疏水亲水参数导入-----

setwd("~/misinteration")
hyb <- read.csv("AA_hyb.csv",header = F, stringsAsFactors = F)
hyb <- hyb[,-1]

ura3_ACC <- read.csv('ura3_pos_ACC.csv',header = F, stringsAsFactors = F)
names(ura3_ACC) <-c('AA_pos','AA_bef','ACC')
Surface_area <- read.csv("surface_area_AA.csv",header = T)
names(Surface_area) <- c("AA_bef","sur_area")
ura3_RSA<- merge(ura3_ACC,Surface_area,by='AA_bef',all=T)
#ura3_RSA$ACC <-as.numeric(ura3_RSA$ACC)
ura3_RSA$RSA <- ura3_RSA$ACC/ura3_RSA$sur_area

names(hyb) <- c("AA_bef","hyb_bef")
ura3_AA_bef_hyb <- merge(ura3_RSA,hyb,by='AA_bef',all = T)

ura3_AA_mutafter <- read.csv('ura3_AA_after.csv',header = F)
names(ura3_AA_mutafter) <-c('mut','AA_pos','AA_aft')

names(hyb) <- c("AA_aft","hyb_aft")
ura3_AA_mutafter_hyb<- merge(ura3_AA_mutafter,hyb,by='AA_aft',all = T)

ura3_AAba_RSA <- merge(ura3_AA_bef_hyb,ura3_AA_mutafter_hyb,by='AA_pos',all = T)
ura3_AAba_RSA <- na.omit(ura3_AAba_RSA)
ura3_AAba_RSA$hyb_t <- ura3_AAba_RSA$hyb_aft - ura3_AAba_RSA$hyb_bef

# TUY ------

TUY_day3$mut <- as.character(TUY_day3$mut)

TUY_day3$AA_pos <- 
  TUY_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

TUY_day3_RSA <- merge(TUY_day3,ura3_AAba_RSA,by=c('AA_pos','mut'),all = T)
TUY_day3_RSA <- na.omit(TUY_day3_RSA)

# misinteraction

TUY_day3_geno <- TUY_day3_RSA[,c(1,2,3,8,12)]

TUY_day3_surfaced <- filter(TUY_day3_geno, RSA > 0.2)
TUY_day3_surfaced <- TUY_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TUY=mean(fitness_TUY))

ggplot(TUY_day3_surfaced,aes(TUY_day3_surfaced$hyb_t,TUY_day3_surfaced$fitness_TUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TUY diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(TUY_day3_surfaced$hyb_t,TUY_day3_surfaced$fitness_TUY,method = "s")

TUY_day3_surfaced$hyb_RSA <- TUY_day3_surfaced$hyb_t * TUY_day3_surfaced$RSA
ggplot(TUY_day3_surfaced,aes(TUY_day3_surfaced$hyb_RSA,TUY_day3_surfaced$fitness_TUY)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TUY_day3_surfaced$hyb_RSA,TUY_day3_surfaced$fitness_TUY,method = "s")


# misfolding 

mut_misfolding_ura3 <- read.csv('mut_misfolding_ura3.csv',header = F)
names(mut_misfolding_ura3) <- c('mut','misfolding')

TUY_day3_misfolding <- merge(TUY_day3_geno,mut_misfolding_ura3,by='mut',all = T)
TUY_day3_misfolding <- na.omit(TUY_day3_misfolding)
TUY_day3_buried_misfolding <- filter(TUY_day3_misfolding, RSA<0.2)

ggplot(TUY_day3_buried_misfolding,aes(TUY_day3_buried_misfolding$misfolding,TUY_day3_buried_misfolding$fitness_TUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(TUY_day3_buried_misfolding$misfolding,TUY_day3_buried_misfolding$fitness_TUY,method = "s")

# IDR 

mut_disorder_ura3 <- read.csv('mut_disorder_ura3_10web.csv',header = F)
names(mut_disorder_ura3) <- c('mut','disorder_length')

TUY_day3_sur_diso <- merge(TUY_day3_geno,mut_disorder_ura3,by='mut',all = T)
TUY_day3_sur_diso <- na.omit(TUY_day3_sur_diso)

ggplot(TUY_day3_sur_diso,aes(TUY_day3_sur_diso$disorder_length,TUY_day3_sur_diso$fitness_TUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TUY_day3_sur_diso$disorder_length,TUY_day3_sur_diso$fitness_TUY,method = "s")

TUY_day3_sur_mean_disor <- TUY_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_TUY = mean(fitness_TUY))

ggplot(TUY_day3_sur_mean_disor,aes(TUY_day3_sur_mean_disor$disorder_length,TUY_day3_sur_mean_disor$fitness_TUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 20, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TUY IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TUY_day3_sur_mean_disor$disorder_length,TUY_day3_sur_mean_disor$fitness_TUY,method = "s")


# AUY ------
AUY_day3$mut <- as.character(AUY_day3$mut)

AUY_day3$AA_pos <- 
  AUY_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

AUY_day3_RSA <- merge(AUY_day3,ura3_AAba_RSA,by=c('AA_pos','mut'),all = T)
AUY_day3_RSA <- na.omit(AUY_day3_RSA)

# misinteraction

AUY_day3_geno <- AUY_day3_RSA[,c(1,2,3,8,12)]

AUY_day3_surfaced <- filter(AUY_day3_geno, RSA > 0.2)
AUY_day3_surfaced <- AUY_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AUY=mean(fitness_AUY))

ggplot(AUY_day3_surfaced,aes(AUY_day3_surfaced$hyb_t,AUY_day3_surfaced$fitness_AUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "AUY diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(AUY_day3_surfaced$hyb_t,AUY_day3_surfaced$fitness_AUY,method = "s")

AUY_day3_surfaced$hyb_RSA <- AUY_day3_surfaced$hyb_t * AUY_day3_surfaced$RSA
ggplot(AUY_day3_surfaced,aes(AUY_day3_surfaced$hyb_RSA,AUY_day3_surfaced$fitness_AUY)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(AUY_day3_surfaced$hyb_RSA,AUY_day3_surfaced$fitness_AUY,method = "s")


# misfolding 

mut_misfolding_ura3 <- read.csv('mut_misfolding_ura3.csv',header = F)
names(mut_misfolding_ura3) <- c('mut','misfolding')

AUY_day3_misfolding <- merge(AUY_day3_geno,mut_misfolding_ura3,by='mut',all = T)
AUY_day3_misfolding <- na.omit(AUY_day3_misfolding)
AUY_day3_buried_misfolding <- filter(AUY_day3_misfolding, RSA<0.2)

ggplot(AUY_day3_buried_misfolding,aes(AUY_day3_buried_misfolding$misfolding,AUY_day3_buried_misfolding$fitness_AUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(AUY_day3_buried_misfolding$misfolding,AUY_day3_buried_misfolding$fitness_AUY,method = "s")

# IDR 

mut_disorder_ura3 <- read.csv('mut_disorder_ura3_10web.csv',header = F)
names(mut_disorder_ura3) <- c('mut','disorder_length')

AUY_day3_sur_diso <- merge(AUY_day3_geno,mut_disorder_ura3,by='mut',all = T)
AUY_day3_sur_diso <- na.omit(AUY_day3_sur_diso)

ggplot(AUY_day3_sur_diso,aes(AUY_day3_sur_diso$disorder_length,AUY_day3_sur_diso$fitness_AUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(AUY_day3_sur_diso$disorder_length,AUY_day3_sur_diso$fitness_AUY,method = "s")

AUY_day3_sur_mean_disor <- AUY_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_AUY = mean(fitness_AUY))

ggplot(AUY_day3_sur_mean_disor,aes(AUY_day3_sur_mean_disor$disorder_length,AUY_day3_sur_mean_disor$fitness_AUY)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 20, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "AUY IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(AUY_day3_sur_mean_disor$disorder_length,AUY_day3_sur_mean_disor$fitness_AUY,method = "s")


# TUS ------
TUS_day3$mut <- as.character(TUS_day3$mut)

TUS_day3$AA_pos <- 
  TUS_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

TUS_day3_RSA <- merge(TUS_day3,ura3_AAba_RSA,by=c('AA_pos','mut'),all = T)
TUS_day3_RSA <- na.omit(TUS_day3_RSA)

# misinteraction

TUS_day3_geno <- TUS_day3_RSA[,c(1,2,3,8,12)]

TUS_day3_surfaced <- filter(TUS_day3_geno, RSA > 0.2)
TUS_day3_surfaced <- TUS_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TUS=mean(fitness_TUS))

ggplot(TUS_day3_surfaced,aes(TUS_day3_surfaced$hyb_t,TUS_day3_surfaced$fitness_TUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TUS diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(TUS_day3_surfaced$hyb_t,TUS_day3_surfaced$fitness_TUS,method = "s")

TUS_day3_surfaced$hyb_RSA <- TUS_day3_surfaced$hyb_t * TUS_day3_surfaced$RSA
ggplot(TUS_day3_surfaced,aes(TUS_day3_surfaced$hyb_RSA,TUS_day3_surfaced$fitness_TUS)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TUS_day3_surfaced$hyb_RSA,TUS_day3_surfaced$fitness_TUS,method = "s")


# misfolding 

mut_misfolding_ura3 <- read.csv('mut_misfolding_ura3.csv',header = F)
names(mut_misfolding_ura3) <- c('mut','misfolding')

TUS_day3_misfolding <- merge(TUS_day3_geno,mut_misfolding_ura3,by='mut',all = T)
TUS_day3_misfolding <- na.omit(TUS_day3_misfolding)
TUS_day3_buried_misfolding <- filter(TUS_day3_misfolding, RSA<0.2)

ggplot(TUS_day3_buried_misfolding,aes(TUS_day3_buried_misfolding$misfolding,TUS_day3_buried_misfolding$fitness_TUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(TUS_day3_buried_misfolding$misfolding,TUS_day3_buried_misfolding$fitness_TUS,method = "s")

# IDR 

mut_disorder_ura3 <- read.csv('mut_disorder_ura3_10web.csv',header = F)
names(mut_disorder_ura3) <- c('mut','disorder_length')

TUS_day3_sur_diso <- merge(TUS_day3_geno,mut_disorder_ura3,by='mut',all = T)
TUS_day3_sur_diso <- na.omit(TUS_day3_sur_diso)

ggplot(TUS_day3_sur_diso,aes(TUS_day3_sur_diso$disorder_length,TUS_day3_sur_diso$fitness_TUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TUS_day3_sur_diso$disorder_length,TUS_day3_sur_diso$fitness_TUS,method = "s")

TUS_day3_sur_mean_disor <- TUS_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_TUS = mean(fitness_TUS))

ggplot(TUS_day3_sur_mean_disor,aes(TUS_day3_sur_mean_disor$disorder_length,TUS_day3_sur_mean_disor$fitness_TUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 20, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TUS IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TUS_day3_sur_mean_disor$disorder_length,TUS_day3_sur_mean_disor$fitness_TUS,method = "s")



# AUS ------
AUS_day3$mut <- as.character(AUS_day3$mut)

AUS_day3$AA_pos <- 
  AUS_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

AUS_day3_RSA <- merge(AUS_day3,ura3_AAba_RSA,by=c('AA_pos','mut'),all = T)
AUS_day3_RSA <- na.omit(AUS_day3_RSA)

# misinteraction

AUS_day3_geno <- AUS_day3_RSA[,c(1,2,3,8,12)]

AUS_day3_surfaced <- filter(AUS_day3_geno, RSA > 0.2)
AUS_day3_surfaced <- AUS_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AUS=mean(fitness_AUS))

ggplot(AUS_day3_surfaced,aes(AUS_day3_surfaced$hyb_t,AUS_day3_surfaced$fitness_AUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "AUS diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(AUS_day3_surfaced$hyb_t,AUS_day3_surfaced$fitness_AUS,method = "s")

AUS_day3_surfaced$hyb_RSA <- AUS_day3_surfaced$hyb_t * AUS_day3_surfaced$RSA
ggplot(AUS_day3_surfaced,aes(AUS_day3_surfaced$hyb_RSA,AUS_day3_surfaced$fitness_AUS)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(AUS_day3_surfaced$hyb_RSA,AUS_day3_surfaced$fitness_AUS,method = "s")


# misfolding 

mut_misfolding_ura3 <- read.csv('mut_misfolding_ura3.csv',header = F)
names(mut_misfolding_ura3) <- c('mut','misfolding')

AUS_day3_misfolding <- merge(AUS_day3_geno,mut_misfolding_ura3,by='mut',all = T)
AUS_day3_misfolding <- na.omit(AUS_day3_misfolding)
AUS_day3_buried_misfolding <- filter(AUS_day3_misfolding, RSA<0.2)

ggplot(AUS_day3_buried_misfolding,aes(AUS_day3_buried_misfolding$misfolding,AUS_day3_buried_misfolding$fitness_AUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(AUS_day3_buried_misfolding$misfolding,AUS_day3_buried_misfolding$fitness_AUS,method = "s")

# IDR 

mut_disorder_ura3 <- read.csv('mut_disorder_ura3_10web.csv',header = F)
names(mut_disorder_ura3) <- c('mut','disorder_length')

AUS_day3_sur_diso <- merge(AUS_day3_geno,mut_disorder_ura3,by='mut',all = T)
AUS_day3_sur_diso <- na.omit(AUS_day3_sur_diso)

ggplot(AUS_day3_sur_diso,aes(AUS_day3_sur_diso$disorder_length,AUS_day3_sur_diso$fitness_AUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(AUS_day3_sur_diso$disorder_length,AUS_day3_sur_diso$fitness_AUS,method = "s")

AUS_day3_sur_mean_disor <- AUS_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_AUS = mean(fitness_AUS))

ggplot(AUS_day3_sur_mean_disor,aes(AUS_day3_sur_mean_disor$disorder_length,AUS_day3_sur_mean_disor$fitness_AUS)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 20, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "AUS IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(AUS_day3_sur_mean_disor$disorder_length,AUS_day3_sur_mean_disor$fitness_AUS,method = "s")


# TUU ------
TUU_day3$mut <- as.character(TUU_day3$mut)

TUU_day3$AA_pos <- 
  TUU_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

TUU_day3_RSA <- merge(TUU_day3,ura3_AAba_RSA,by=c('AA_pos','mut'),all = T)
TUU_day3_RSA <- na.omit(TUU_day3_RSA)

# misinteraction

TUU_day3_geno <- TUU_day3_RSA[,c(1,2,3,8,12)]

TUU_day3_surfaced <- filter(TUU_day3_geno, RSA > 0.2)
TUU_day3_surfaced <- TUU_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TUU=mean(fitness_TUU))

ggplot(TUU_day3_surfaced,aes(TUU_day3_surfaced$hyb_t,TUU_day3_surfaced$fitness_TUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TUU diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(TUU_day3_surfaced$hyb_t,TUU_day3_surfaced$fitness_TUU,method = "s")

TUU_day3_surfaced$hyb_RSA <- TUU_day3_surfaced$hyb_t * TUU_day3_surfaced$RSA
ggplot(TUU_day3_surfaced,aes(TUU_day3_surfaced$hyb_RSA,TUU_day3_surfaced$fitness_TUU)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TUU_day3_surfaced$hyb_RSA,TUU_day3_surfaced$fitness_TUU,method = "s")


# misfolding 

mut_misfolding_ura3 <- read.csv('mut_misfolding_ura3.csv',header = F)
names(mut_misfolding_ura3) <- c('mut','misfolding')

TUU_day3_misfolding <- merge(TUU_day3_geno,mut_misfolding_ura3,by='mut',all = T)
TUU_day3_misfolding <- na.omit(TUU_day3_misfolding)
TUU_day3_buried_misfolding <- filter(TUU_day3_misfolding, RSA<0.2)

ggplot(TUU_day3_buried_misfolding,aes(TUU_day3_buried_misfolding$misfolding,TUU_day3_buried_misfolding$fitness_TUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(TUU_day3_buried_misfolding$misfolding,TUU_day3_buried_misfolding$fitness_TUU,method = "s")

# IDR 

mut_disorder_ura3 <- read.csv('mut_disorder_ura3_10web.csv',header = F)
names(mut_disorder_ura3) <- c('mut','disorder_length')

TUU_day3_sur_diso <- merge(TUU_day3_geno,mut_disorder_ura3,by='mut',all = T)
TUU_day3_sur_diso <- na.omit(TUU_day3_sur_diso)

ggplot(TUU_day3_sur_diso,aes(TUU_day3_sur_diso$disorder_length,TUU_day3_sur_diso$fitness_TUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TUU_day3_sur_diso$disorder_length,TUU_day3_sur_diso$fitness_TUU,method = "s")

TUU_day3_sur_mean_disor <- TUU_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness_TUU = mean(fitness_TUU))

ggplot(TUU_day3_sur_mean_disor,aes(TUU_day3_sur_mean_disor$disorder_length,TUU_day3_sur_mean_disor$fitness_TUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 20, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TUU IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TUU_day3_sur_mean_disor$disorder_length,TUU_day3_sur_mean_disor$fitness_TUU,method = "s")


# AUU ------
AUU_day3$mut <- as.character(AUU_day3$mut)

AUU_day3$AA_pos <- 
  AUU_day3$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){ceiling(as.numeric(paste(x[2:4],collapse="",sep=""))/3)}else if(length(x)==4)
  {ceiling(as.numeric(paste(x[2:3],collapse="",sep=""))/3)}else if(length(x)==3){ceiling(as.numeric(x[2])/3)}}) %>% 
  unlist()

AUU_day3_RSA <- merge(AUU_day3,ura3_AAba_RSA,by=c('AA_pos','mut'),all = T)
AUU_day3_RSA <- na.omit(AUU_day3_RSA)

# misinteraction

AUU_day3_geno <- AUU_day3_RSA[,c(1,2,3,8,12)]

AUU_day3_surfaced <- filter(AUU_day3_geno, RSA > 0.2)
AUU_day3_surfaced <- AUU_day3_surfaced %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AUU=mean(fitness_AUU))

ggplot(AUU_day3_surfaced,aes(AUU_day3_surfaced$hyb_t,AUU_day3_surfaced$fitness_AUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "AUU diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))
cor.test(AUU_day3_surfaced$hyb_t,AUU_day3_surfaced$fitness_AUU,method = "s")

AUU_day3_surfaced$hyb_RSA <- AUU_day3_surfaced$hyb_t * AUU_day3_surfaced$RSA
ggplot(AUU_day3_surfaced,aes(AUU_day3_surfaced$hyb_RSA,AUU_day3_surfaced$fitness_AUU)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(AUU_day3_surfaced$hyb_RSA,AUU_day3_surfaced$fitness_AUU,method = "s")


# misfolding 

mut_misfolding_ura3 <- read.csv('mut_misfolding_ura3.csv',header = F)
names(mut_misfolding_ura3) <- c('mut','misfolding')

AUU_day3_misfolding <- merge(AUU_day3_geno,mut_misfolding_ura3,by='mut',all = T)
AUU_day3_misfolding <- na.omit(AUU_day3_misfolding)
AUU_day3_buried_misfolding <- filter(AUU_day3_misfolding, RSA<0.2)

ggplot(AUU_day3_buried_misfolding,aes(AUU_day3_buried_misfolding$misfolding,AUU_day3_buried_misfolding$fitness_AUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("misfolding") + ylab("fitness")+ style.print()

cor.test(AUU_day3_buried_misfolding$misfolding,AUU_day3_buried_misfolding$fitness_AUU,method = "s")

# IDR 

mut_disorder_ura3 <- read.csv('mut_disorder_ura3_10web.csv',header = F)
names(mut_disorder_ura3) <- c('mut','disorder_length')

AUU_day3_sur_diso <- merge(AUU_day3_geno,mut_disorder_ura3,by='mut',all = T)
AUU_day3_sur_diso <- na.omit(AUU_day3_sur_diso)

ggplot(AUU_day3_sur_diso,aes(AUU_day3_sur_diso$disorder_length,AUU_day3_sur_diso$fitness_AUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(AUU_day3_sur_diso$disorder_length,AUU_day3_sur_diso$fitness_AUU,method = "s")

AUU_day3_sur_mean_disor <- AUU_day3_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
dplyr::summarise(fitness_AUU = mean(fitness_AUU))

ggplot(AUU_day3_sur_mean_disor,aes(AUU_day3_sur_mean_disor$disorder_length,AUU_day3_sur_mean_disor$fitness_AUU)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P=   rho= ', x = 20, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "AUU IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(AUU_day3_sur_mean_disor$disorder_length,AUU_day3_sur_mean_disor$fitness_AUU,method = "s")






