TGY_day0wt <- read.csv("TGY_day0_wt.csv",header = F)
names(TGY_day0wt) <- c("umi","day0")

TGY_day7wt <- read.csv("TGY_day7_wt.csv",header = F)
names(TGY_day7wt) <- c("umi","day7")

TGY_day0single <- read.csv("TGY_day0_single.csv",header = F)
names(TGY_day0single) <- c("umi","pos","mut","DNA","AA_pos","mut_AA","day0")

TGY_day7single <- read.csv("TGY_day7_single.csv",header = F)
names(TGY_day7single) <- c("umi","pos","mut","DNA","AA_pos","mut_AA","day7")

TGY_day07_wt <- merge(TGY_day0wt,TGY_day7wt,by="umi",all = T)
TGY_day07_wt[is.na(TGY_day07_wt)]<-0
TGY_day07_wt <- TGY_day07_wt %>%
  filter(day0>0) %>%
  mutate(wt_ratio = day7/day0)

quantile(TGY_day07_wt$wt_ratio,0.975)
quantile(TGY_day07_wt$wt_ratio,0.025)
TGY_day07_wt <- filter(TGY_day07_wt,wt_ratio<quantile(TGY_day07_wt$wt_ratio,0.975) & wt_ratio>quantile(TGY_day07_wt$wt_ratio,0.025))

sum(TGY_day07_wt$day0)
18585670
sum(TGY_day07_wt$day7)
27886116

TGY_day07_single <- merge(TGY_day0single,TGY_day7single,by=c("umi","pos","mut","DNA","AA_pos","mut_AA"),all = T)
TGY_day07_single[is.na(TGY_day07_single)]<-0

SC-74.7

# 氨基酸顺序和突变前AA疏水性-------
aa_que <- read.csv("AA_que.csv",header = F,stringsAsFactors = F,colClasses = "character")
aa_que <- t(aa_que)
aa_que <- as.data.frame(aa_que)
names(aa_que) <- c("AA_pos","AA_bef")   #没问题 aa_que是225个氨基酸顺序

aa_que$AA_pos <- as.factor(aa_que$AA_pos)%>%as.character()
aa_que$AA_bef <- as.factor(aa_que$AA_bef)%>%as.character()

hyb2 <- read.csv("AA_hyb.csv",header = F, stringsAsFactors = F)
names(hyb2) <- c("AA","AA_bef","hydrophobicity")

hyb2$AA_bef<-as.character(hyb2$AA_bef)

aa_que1 <- merge(aa_que,hyb2,by = "AA_bef",all=T)

aa_que1$AA_pos <- as.numeric(aa_que1$AA_pos)

#引入疏水性表------
hyb <- read.csv("AA_hyb.csv",header = F, stringsAsFactors = F)
names(hyb) <- c("AA","mut_AA","hydrophobicity")
TGY_day07_single <- merge(TGY_day07_single,hyb,by = "mut_AA",all=T)


# 删去顺序里没有的并merge------ (后面gfp做了更改)
w<-which(TGY_day07_single$AA_pos=="1" | TGY_day07_single$AA_pos=="65" | TGY_day07_single$AA_pos=="66" | TGY_day07_single$AA_pos=="67" |TGY_day07_single$AA_pos>=230 | TGY_day07_single$mut_AA =="*")
TGY_day07_single1 <- TGY_day07_single[-w,]

TGY_day07_single1$AA_pos <- ifelse(TGY_day07_single1$AA_pos <65, TGY_day07_single1$AA_pos - 1, TGY_day07_single1$AA_pos-4)

TGY07_ba_hyb <- merge(TGY_day07_single1,aa_que1,by="AA_pos",all = T)
TGY07_ba_hyb$hydrophobicity.x <- as.numeric(TGY07_ba_hyb$hydrophobicity.x)
TGY07_ba_hyb$hydrophobicity.y <- as.numeric(TGY07_ba_hyb$hydrophobicity.y)
TGY07_ba_hyb$hyb_t <- TGY07_ba_hyb$hydrophobicity.x - TGY07_ba_hyb$hydrophobicity.y

# 引入ACC/RSA-----

ACC <- read.csv("ACC.csv",header = T, stringsAsFactors = F)
names(ACC) <- c("AA_bef", "ACC")
ACC$AA_pos <- 1:225

RSA_pre<- read.csv("ACC.csv",header = T)
RSA<-as.data.frame(RSA_pre)
RSA$Num<-1:225
names(RSA)<-c("AA","ACC","AA_pos")
Surface_area <- read.csv("surface_area_AA.csv",header = T)
names(Surface_area) <- c("AA","sur_area")
RSA$AA <-as.character(RSA$AA)
Surface_area$AA <- as.character(Surface_area$AA)
RSA1 <- merge(RSA,Surface_area,by="AA")
RSA1$RSA <- RSA1$ACC/RSA1$sur_area

RSA2 <- RSA1[,-c(2,4)]
names(RSA2) <-c("AA_bef","AA_pos","RSA")

# ACC/RSA merge进大数据集------

TGY07_hyb_RA <- merge(TGY07_ba_hyb,ACC,by=c("AA_bef","AA_pos"),all = T)
TGY07_hyb_RA <- merge(TGY07_hyb_RA,RSA2,by=c("AA_bef","AA_pos"),all = T)

TGY_type <-TGY07_hyb_RA
TGY_type$type <- ifelse(TGY_type$RSA>=0.2, "surfaced", "buried")

TGY_type %>% 
  ggplot(aes(x = fitness,color=type))+
  geom_density() + style.print()


#200531 突变geno一样的合并-----

TGY_geno <- TGY07_hyb_RA %>%
  dplyr::group_by(AA_bef,AA_pos,mut_AA,pos,mut,DNA,AA.x,AA.y,hydrophobicity.x,hydrophobicity.y,hyb_t,ACC,RSA) %>%
  dplyr::summarise(day0 = sum(day0),day7 = sum(day7))

TGY_geno <- filter(TGY_geno,day0 >100)

TGY_geno <- TGY_geno %>%
  mutate(fitness = (day7 / day0 / sum(TGY_day07_wt$day7) * sum(TGY_day07_wt$day0)) ** (1/106.1))


#RSA0.2------
TGY_RSA0.2_surfaced <- filter(TGY_geno, RSA > 0.2 & hyb_t!=0)

TGY_RSA0.2_buried <- filter(TGY_geno, RSA < 0.2 & hyb_t != 0)

wilcox.test(TGY_RSA0.2_buried$fitness, TGY_RSA0.2_surfaced$fitness)
t.test(TGY_RSA0.2_buried$fitness, TGY_RSA0.2_surfaced$fitness)

TGY_RSA0.2_surfaced$hyb_ACC <- TGY_RSA0.2_surfaced$hyb_t * TGY_RSA0.2_surfaced$ACC
ggplot(TGY_RSA0.2_surfaced,aes(TGY_RSA0.2_surfaced$hyb_ACC,TGY_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*ACC") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*ACC

cor.test(TGY_RSA0.2_surfaced$hyb_ACC,TGY_RSA0.2_surfaced$fitness,method = "s")

ggplot(TGY_RSA0.2_surfaced,aes(TGY_RSA0.2_surfaced$hyb_t,TGY_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0014  rho= -0.12', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TGY diff_hyb sufarced")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))

cor.test(TGY_RSA0.2_surfaced$hyb_t,TGY_RSA0.2_surfaced$fitness,method = "s")

TGY_RSA0.2_surfaced$hyb_RSA <- TGY_RSA0.2_surfaced$hyb_t * TGY_RSA0.2_surfaced$RSA
ggplot(TGY_RSA0.2_surfaced,aes(TGY_RSA0.2_surfaced$hyb_RSA,TGY_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TGY_RSA0.2_surfaced$hyb_RSA,TGY_RSA0.2_surfaced$fitness,method = "s")

# buried 0.2------

ggplot(TGY_RSA0.2_buried,aes(TGY_RSA0.2_buried$hyb_t,TGY_RSA0.2_buried$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1427  rho= 0.05', x = 5, y = 1.01, size=7, colour ='#000000')+
  labs(title = "TGY diff_hyb buried")+ 
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, face = "bold"))

cor.test(TGY_RSA0.2_buried$hyb_t,TGY_RSA0.2_buried$fitness,method = "s")

#只看hyb>0的-----

TGY_0.2_surfaced_hy_more0 <- filter(TGY_geno, RSA > 0.2 & hyb_t>0)

ggplot(TGY_0.2_surfaced_hy_more0,aes(TGY_0.2_surfaced_hy_more0$hyb_t,TGY_0.2_surfaced_hy_more0$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGY_0.2_surfaced_hy_more0$hyb_t,TGY_0.2_surfaced_hy_more0$fitness,method = "s")

#把小于0的取绝对值，画出差距图------

hyb_RSA0.2 <- filter(TGY_geno, RSA > 0.2 & hyb_t>0)
hyl_RSA0.2 <- filter(TGY_geno, RSA > 0.2 & hyb_t<0)

hyl_RSA0.2_abs <- hyl_RSA0.2
hyl_RSA0.2_abs$hyb_t <- hyl_RSA0.2_abs$hyb_t * -1

hyb_RSA0.2$type <- "hyb"
hyl_RSA0.2_abs$type <- "hyl"

hyb_hyl_abs <- rbind(hyb_RSA0.2, hyl_RSA0.2_abs)

ggplot(hyb_hyl_abs,aes(hyb_hyl_abs$hyb_t,hyb_hyl_abs$fitness, fill= type)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness") + style.print()   #绝对值线

wilcox.test(hyb_RSA0.2$fitness, hyl_RSA0.2_abs$fitness)

# 上图随着x轴增加，两条线fitness中值和平均值的差距变化-----

data <- lapply(c(1:8),function(x){
  cutoff2_hyb <- filter(hyb_hyl_abs, hyb_t >x & type == "hyb")
  cutoff2_hyl <- filter(hyb_hyl_abs, hyb_t >x & type == "hyl")
  data1 <- data.frame(cutoff = x, hyb_mean = mean(cutoff2_hyb$fitness), hyl_mean = mean(cutoff2_hyl$fitness), hyb_med = median(cutoff2_hyb$fitness), hyl_med = median(cutoff2_hyl$fitness))
}) %>% rbind.fill()

data2 <- melt(data, id.vars = c("cutoff"))

data2_mean <- filter(data2, variable == "hyl_mean" |  variable =="hyb_mean")

ggplot(data2_mean, aes(cutoff,value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity") + coord_cartesian(ylim = c(1, 1.005)) + 
  theme(legend.position='none')+ xlab("changed hydrophobicity cutoff") + ylab("fitness mean") + style.print() 

data2_med <- filter(data2, variable == "hyl_med" |  variable =="hyb_med")

ggplot(data2_med, aes(cutoff,value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity") + coord_cartesian(ylim = c(1, 1.005)) + 
  theme(legend.position='none') + xlab("changed hydrophobicity cutoff") + ylab("fitness med") + style.print()

#把两端的值去掉（-4到4）------

TGY_RSA0.2_cutsides4 <- filter(TGY_geno, RSA > 0.2 & hyb_t>-4 & hyb_t<4 & hyb_t!=0)

ggplot(TGY_RSA0.2_cutsides4,aes(TGY_RSA0.2_cutsides4$hyb_t,TGY_RSA0.2_cutsides4$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGY_RSA0.2_cutsides4$hyb_t,TGY_RSA0.2_cutsides4$fitness,method = "s")

#下列为重复（新数据除两端）
#————————————————————————————

hyb_RSA0.2 <- filter(TGY_geno, RSA > 0.2 & hyb_t<4 & hyb_t>0)
hyl_RSA0.2 <- filter(TGY_geno, RSA > 0.2 & hyb_t>-4 & hyb_t<0)

hyl_RSA0.2_abs <- hyl_RSA0.2
hyl_RSA0.2_abs$hyb_t <- hyl_RSA0.2_abs$hyb_t * -1

hyb_RSA0.2$type <- "hyb"
hyl_RSA0.2_abs$type <- "hyl"

hyb_hyl_abs <- rbind(hyb_RSA0.2, hyl_RSA0.2_abs)

ggplot(hyb_hyl_abs,aes(hyb_hyl_abs$hyb_t,hyb_hyl_abs$fitness, fill= type)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness") + style.print()   #绝对值线

t.test(hyb_RSA0.2$fitness, hyl_RSA0.2_abs$fitness)

data <- lapply(c(1:4),function(x){
  cutoff2_hyb <- filter(hyb_hyl_abs, hyb_t >x & type == "hyb")
  cutoff2_hyl <- filter(hyb_hyl_abs, hyb_t >x & type == "hyl")
  data1 <- data.frame(cutoff = x, hyb_mean = mean(cutoff2_hyb$fitness), hyl_mean = mean(cutoff2_hyl$fitness), hyb_med = median(cutoff2_hyb$fitness), hyl_med = median(cutoff2_hyl$fitness))
}) %>% rbind.fill()

data2 <- melt(data, id.vars = c("cutoff"))

data2_mean <- filter(data2, variable == "hyl_mean" |  variable =="hyb_mean")

ggplot(data2_mean, aes(cutoff,value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity") + coord_cartesian(ylim = c(1, 1.005)) + 
  theme(legend.position='none')+ xlab("changed hydrophobicity cutoff") + ylab("fitness mean") + style.print() 

data2_med <- filter(data2, variable == "hyl_med" |  variable =="hyb_med")

ggplot(data2_med, aes(cutoff,value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity") + coord_cartesian(ylim = c(1, 1.005)) + 
  theme(legend.position='none') + xlab("changed hydrophobicity cutoff") + ylab("fitness med") + style.print()

#—————————————————————————————————————————————


#看疏水性符号变化了的 选取-7到7------

hyb_changed_RSA0.2 <- filter(TGY_geno, RSA>0.2 & hydrophobicity.x<0 & hydrophobicity.y>0 & hyb_t <7 & hyb_t>-7 |RSA>0.2 & hydrophobicity.x>0 & hydrophobicity.y<0 & hyb_t <7 & hyb_t>-7)

ggplot(hyb_changed_RSA0.2,aes(hyb_changed_RSA0.2$hyb_t,hyb_changed_RSA0.2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(hyb_changed_RSA0.2$hyb_t,hyb_changed_RSA0.2$fitness,method = "s")


#看亲水到疏水+没那么亲水，疏水到没那么疏水+亲水 t.test-----


hyb_changedless_RSA0.2 <- filter(TGY_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y |RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

ggplot(hyb_changedless_RSA0.2,aes(hyb_changedless_RSA0.2$hyb_t,hyb_changedless_RSA0.2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(hyb_changedless_RSA0.2$hyb_t,hyb_changedless_RSA0.2$fitness,method = "s")


aa_RSA0.2_sur_qin_shu <- filter(TGY_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

aa_RSA0.2_sur_shu_qin_lesshyb <- filter(TGY_geno, RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

t.test(aa_RSA0.2_sur_qin_shu$fitness, aa_RSA0.2_sur_shu_qin_lesshyb$fitness)




#200527 把相同位置的蛋白合在一起（fitness求均值），看位置突变与fitness的关系-----

TGY_pos <- TGY_type %>%
  dplyr::group_by(AA_pos,AA_bef,ACC,RSA,type) %>%
  dplyr::summarise(day0 = mean(day0),day7 = mean(day7))

TGY_pos <- TGY_type %>%
  dplyr::group_by(pos,type) %>%
  dplyr::summarise(day0 = mean(day0),day7 = mean(day7))


TGY_pos <- TGY_pos %>%
  mutate(fitness = (day7 / day0 / 6284442 * 8758950) ** (1/106.1))

ggplot(TGY_pos,aes(x=pos,y=fitness,colour=type))+
  geom_point()


#200619 散点图改成分组图-----

TGY_hybt_fitness <- TGY_RSA0.2_surfaced[,c(11,16)]
TGY_hybt_fitness1 <- TGY_hybt_fitness[sort(TGY_hybt_fitness$hyb_t, index.return = T)$ix,]
TGY_hybt_fitness1$seq <- 1:721

TGY_hybt_fitness1[,4] <- "group"

#十组------

for (i in 1:nrow(TGY_hybt_fitness1)){
  if(TGY_hybt_fitness1[i,3]<72){
  TGY_hybt_fitness1[i,4] <-"0-72"
}else if(TGY_hybt_fitness1[i,3]<144){
  TGY_hybt_fitness1[i,4] <-"72-144"
}else if(TGY_hybt_fitness1[i,3]<216){
  TGY_hybt_fitness1[i,4] <-"144-216"
}else if(TGY_hybt_fitness1[i,3]<288){
  TGY_hybt_fitness1[i,4] <-"216-288"
}else if(TGY_hybt_fitness1[i,3]<360){
  TGY_hybt_fitness1[i,4] <-"288-360"
}else if(TGY_hybt_fitness1[i,3]<432){
  TGY_hybt_fitness1[i,4] <-"360-432"
}else if(TGY_hybt_fitness1[i,3]<504){
  TGY_hybt_fitness1[i,4] <-"432-504"
}else if(TGY_hybt_fitness1[i,3]<576){
  TGY_hybt_fitness1[i,4] <-"504-576"
}else if(TGY_hybt_fitness1[i,3]<648){
  TGY_hybt_fitness1[i,4] <-"576-648"
}else{
  TGY_hybt_fitness1[i,4] <-"648-721"
}}

names(TGY_hybt_fitness1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness1$group<- factor(TGY_hybt_fitness1$group,levels=c("0-72","72-144","144-216","216-288","288-360","360-432","432-504","504-576","576-648","648-721"))

TGY_hybt_fitness1 %>%
 ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) 


#五组------

for (i in 1:nrow(TGY_hybt_fitness1)){
  if(TGY_hybt_fitness1[i,3]<144){
    TGY_hybt_fitness1[i,4] <-"0-144"
  }else if(TGY_hybt_fitness1[i,3]<288){
    TGY_hybt_fitness1[i,4] <-"144-288"
  }else if(TGY_hybt_fitness1[i,3]<432){
    TGY_hybt_fitness1[i,4] <-"288-432"
  }else if(TGY_hybt_fitness1[i,3]<576){
    TGY_hybt_fitness1[i,4] <-"432-576"
  }else{
    TGY_hybt_fitness1[i,4] <-"576-721"
  }}

names(TGY_hybt_fitness1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness1$group<- factor(TGY_hybt_fitness1$group,levels=c("0-144","144-288","288-432","432-576","576-721"))

means <- aggregate(fitness ~  group, TGY_hybt_fitness1, mean)
medians <- aggregate(fitness ~  group, TGY_hybt_fitness1, median) %>%
  mutate(size = round(1-fitness,6))

TGY_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() + ylim(0.99,1.01)+
  stat_summary(fun.y=median, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = medians, aes(label = size, y = fitness + 0.01))+
  labs(title = "TGY Hyb_t",subtitle = "Median 5groups,  down 317%")

#四组-------

for (i in 1:nrow(TGY_hybt_fitness1)){
  if(TGY_hybt_fitness1[i,3]<180){
    TGY_hybt_fitness1[i,4] <-"0-180"
  }else if(TGY_hybt_fitness1[i,3]<360){
    TGY_hybt_fitness1[i,4] <-"180-360"
  }else if(TGY_hybt_fitness1[i,3]<540){
    TGY_hybt_fitness1[i,4] <-"360-540"
  }else{
    TGY_hybt_fitness1[i,4] <-"540-721"
  }}

names(TGY_hybt_fitness1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness1, median)%>%
  mutate(size = round(1-fitness,6))
means <- aggregate(fitness ~  group, TGY_hybt_fitness1, mean)%>%
  mutate(size = round(1-fitness,6))

TGY_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=median, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = medians, aes(label = size, y = fitness + 0.002)) + ylim(0.99,1.01)+
  labs(title = "TGY Hyb_t",subtitle = "Median 4groups,  down 127%")

(medians[1,2]-medians[4,2])/medians[1,2]

#(插入)判断正态性-----
shapiro.test(TGY_hybt_fitness1$fitness)
ks.test(TGY_hybt_fitness1$fitness, "pnorm", mean = mean(TGY_hybt_fitness1$fitness), sd =  sqrt(var(TGY_hybt_fitness1$fitness)))

#三组------

for (i in 1:nrow(TGY_hybt_fitness1)){
  if(TGY_hybt_fitness1[i,3]<240){
    TGY_hybt_fitness1[i,4] <-"0-240"
  }else if(TGY_hybt_fitness1[i,3]<480){
    TGY_hybt_fitness1[i,4] <-"240-480"
  }else{
    TGY_hybt_fitness1[i,4] <-"480-721"
  }}

names(TGY_hybt_fitness1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = size, y = fitness + 0.002))+ ylim(0.99,1.01) +
  labs(title = "TGY Hyb_t",subtitle = "Mean 3groups,  down 317%") +
  theme(plot.title = element_text(size = 20))

(means[1,2]-means[3,2])/means[1,2]




#看ACC的下降值------


TGY_hybt_fitness_ACC <- TGY_RSA0.2_surfaced[,c(17,16)]
TGY_hybt_fitness_ACC1 <- TGY_hybt_fitness_ACC[sort(TGY_hybt_fitness_ACC$hyb_ACC, index.return = T)$ix,]
TGY_hybt_fitness_ACC1$seq <- 1:721

TGY_hybt_fitness_ACC1[,4] <- "group"

#十组------

for (i in 1:nrow(TGY_hybt_fitness_ACC1)){
  if(TGY_hybt_fitness_ACC1[i,3]<72){
    TGY_hybt_fitness_ACC1[i,4] <-"0-72"
  }else if(TGY_hybt_fitness_ACC1[i,3]<144){
    TGY_hybt_fitness_ACC1[i,4] <-"72-144"
  }else if(TGY_hybt_fitness_ACC1[i,3]<216){
    TGY_hybt_fitness_ACC1[i,4] <-"144-216"
  }else if(TGY_hybt_fitness_ACC1[i,3]<288){
    TGY_hybt_fitness_ACC1[i,4] <-"216-288"
  }else if(TGY_hybt_fitness_ACC1[i,3]<360){
    TGY_hybt_fitness_ACC1[i,4] <-"288-360"
  }else if(TGY_hybt_fitness_ACC1[i,3]<432){
    TGY_hybt_fitness_ACC1[i,4] <-"360-432"
  }else if(TGY_hybt_fitness_ACC1[i,3]<504){
    TGY_hybt_fitness_ACC1[i,4] <-"432-504"
  }else if(TGY_hybt_fitness_ACC1[i,3]<576){
    TGY_hybt_fitness_ACC1[i,4] <-"504-576"
  }else if(TGY_hybt_fitness_ACC1[i,3]<648){
    TGY_hybt_fitness_ACC1[i,4] <-"576-648"
  }else{
    TGY_hybt_fitness_ACC1[i,4] <-"648-721"
  }}

names(TGY_hybt_fitness_ACC1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness_ACC1$group<- factor(TGY_hybt_fitness_ACC1$group,levels=c("0-72","72-144","144-216","216-288","288-360","360-432","432-504","504-576","576-648","648-721"))

TGY_hybt_fitness_ACC1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot()


#五组------

for (i in 1:nrow(TGY_hybt_fitness_ACC1)){
  if(TGY_hybt_fitness_ACC1[i,3]<144){
    TGY_hybt_fitness_ACC1[i,4] <-"0-144"
  }else if(TGY_hybt_fitness_ACC1[i,3]<288){
    TGY_hybt_fitness_ACC1[i,4] <-"144-288"
  }else if(TGY_hybt_fitness_ACC1[i,3]<432){
    TGY_hybt_fitness_ACC1[i,4] <-"288-432"
  }else if(TGY_hybt_fitness_ACC1[i,3]<576){
    TGY_hybt_fitness_ACC1[i,4] <-"432-576"
  }else{
    TGY_hybt_fitness_ACC1[i,4] <-"576-721"
  }}

names(TGY_hybt_fitness_ACC1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness_ACC1$group<- factor(TGY_hybt_fitness_ACC1$group,levels=c("0-144","144-288","288-432","432-576","576-721"))

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness_ACC1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=median, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = medians, aes(label = size, y = fitness + 0.002)) + ylim(0.99,1.008) +
  labs(title = "TGY Hyb_t*ACC",subtitle = "Median 5roups,  down 317%")+
  theme(plot.title = element_text(size = 20))

#用errorbar------

TGY_hybt_fitness_ACC2 <- TGY_hybt_fitness_ACC1 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(fitness = mean(fitness))
sds <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, sd)

TGY_hybt_fitness_ACC2 <- cbind(TGY_hybt_fitness_ACC2,sds)
TGY_hybt_fitness_ACC2 <- TGY_hybt_fitness_ACC2[,c(1,2,4)]
names(TGY_hybt_fitness_ACC2) <- c("group","fitness","sd")

TGY_hybt_fitness_ACC2 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean(fitness)- sd, ymax = mean(fitness)+ sd), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.01))

(means[1,2]-means[5,2])/means[1,2]

#四组-----

for (i in 1:nrow(TGY_hybt_fitness_ACC1)){
  if(TGY_hybt_fitness_ACC1[i,3]<180){
    TGY_hybt_fitness_ACC1[i,4] <-"0-180"
  }else if(TGY_hybt_fitness_ACC1[i,3]<360){
    TGY_hybt_fitness_ACC1[i,4] <-"180-360"
  }else if(TGY_hybt_fitness_ACC1[i,3]<540){
    TGY_hybt_fitness_ACC1[i,4] <-"360-540"
  }else{
    TGY_hybt_fitness_ACC1[i,4] <-"540-721"
  }}

names(TGY_hybt_fitness_ACC1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness_ACC1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = size, y = fitness + 0.002)) + ylim(0.99,1.008) +
  labs(title = "TGY Hyb_t*ACC",subtitle = "Mean 4roups,  down 209%")+
  theme(plot.title = element_text(size = 20))

(means[1,2]-means[4,2])/means[1,2]

#(插入)判断正态性
shapiro.test(TGY_hybt_fitness_ACC1$fitness)
ks.test(TGY_hybt_fitness_ACC1$fitness, "pnorm", mean = mean(TGY_hybt_fitness_ACC1$fitness), sd =  sqrt(var(TGY_hybt_fitness_ACC1$fitness)))

#三组------

for (i in 1:nrow(TGY_hybt_fitness_ACC1)){
  if(TGY_hybt_fitness_ACC1[i,3]<240){
    TGY_hybt_fitness_ACC1[i,4] <-"0-240"
  }else if(TGY_hybt_fitness_ACC1[i,3]<480){
    TGY_hybt_fitness_ACC1[i,4] <-"240-480"
  }else{
    TGY_hybt_fitness_ACC1[i,4] <-"480-721"
  }}

names(TGY_hybt_fitness_ACC1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness_ACC1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = size, y = fitness + 0.002)) +
  ylim(0.99,1.008) +
  labs(title = "TGY Hyb_t*ACC",subtitle = "Mean 3roups,  down 400%")+
  theme(plot.title = element_text(size = 20))

#用errorbar------

TGY_hybt_fitness_ACC2 <- TGY_hybt_fitness_ACC1 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(fitness = mean(fitness))
sds <- aggregate(fitness ~  group, TGY_hybt_fitness_ACC1, sd)

TGY_hybt_fitness_ACC2 <- cbind(TGY_hybt_fitness_ACC2,sds)
TGY_hybt_fitness_ACC2 <- TGY_hybt_fitness_ACC2[,c(1,2,4)]
names(TGY_hybt_fitness_ACC2) <- c("group","fitness","sd")

TGY_hybt_fitness_ACC2 %>%
  ggplot(aes(x=group,y=fitness,group=1)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean(fitness)- sd, ymax = mean(fitness)+ sd), 
              width = 0.2, position = position_dodge(0.9)) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.005))+
  ylim(0.998,1.009) +
  labs(title = "TGY Hyb_t*ACC",subtitle = "Mean 4groups,  down 0.8%")
  

(medians[1,2]-medians[3,2])/medians[1,2]
(means[1,2]-means[3,2])/means[1,2]



#RSA作为系数的图------

TGY_hybt_fitness_RSA <- TGY_RSA0.2_surfaced[,c(18,16)]
TGY_hybt_fitness_RSA1 <- TGY_hybt_fitness_RSA[sort(TGY_hybt_fitness_RSA$hyb_RSA, index.return = T)$ix,]
TGY_hybt_fitness_RSA1$seq <- 1:721

TGY_hybt_fitness_RSA1[,4] <- "group"

#十组------

for (i in 1:nrow(TGY_hybt_fitness_RSA1)){
  if(TGY_hybt_fitness_RSA1[i,3]<72){
    TGY_hybt_fitness_RSA1[i,4] <-"0-72"
  }else if(TGY_hybt_fitness_RSA1[i,3]<144){
    TGY_hybt_fitness_RSA1[i,4] <-"72-144"
  }else if(TGY_hybt_fitness_RSA1[i,3]<216){
    TGY_hybt_fitness_RSA1[i,4] <-"144-216"
  }else if(TGY_hybt_fitness_RSA1[i,3]<288){
    TGY_hybt_fitness_RSA1[i,4] <-"216-288"
  }else if(TGY_hybt_fitness_RSA1[i,3]<360){
    TGY_hybt_fitness_RSA1[i,4] <-"288-360"
  }else if(TGY_hybt_fitness_RSA1[i,3]<432){
    TGY_hybt_fitness_RSA1[i,4] <-"360-432"
  }else if(TGY_hybt_fitness_RSA1[i,3]<504){
    TGY_hybt_fitness_RSA1[i,4] <-"432-504"
  }else if(TGY_hybt_fitness_RSA1[i,3]<576){
    TGY_hybt_fitness_RSA1[i,4] <-"504-576"
  }else if(TGY_hybt_fitness_RSA1[i,3]<648){
    TGY_hybt_fitness_RSA1[i,4] <-"576-648"
  }else{
    TGY_hybt_fitness_RSA1[i,4] <-"648-721"
  }}

names(TGY_hybt_fitness_RSA1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness_RSA1$group<- factor(TGY_hybt_fitness_RSA1$group,levels=c("0-72","72-144","144-216","216-288","288-360","360-432","432-504","504-576","576-648","648-721"))

TGY_hybt_fitness_RSA1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot()


#五组-----

for (i in 1:nrow(TGY_hybt_fitness_RSA1)){
  if(TGY_hybt_fitness_RSA1[i,3]<144){
    TGY_hybt_fitness_RSA1[i,4] <-"0-144"
  }else if(TGY_hybt_fitness_RSA1[i,3]<288){
    TGY_hybt_fitness_RSA1[i,4] <-"144-288"
  }else if(TGY_hybt_fitness_RSA1[i,3]<432){
    TGY_hybt_fitness_RSA1[i,4] <-"288-432"
  }else if(TGY_hybt_fitness_RSA1[i,3]<576){
    TGY_hybt_fitness_RSA1[i,4] <-"432-576"
  }else{
    TGY_hybt_fitness_RSA1[i,4] <-"576-721"
  }}

names(TGY_hybt_fitness_RSA1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness_RSA1$group<- factor(TGY_hybt_fitness_RSA1$group,levels=c("0-144","144-288","288-432","432-576","576-721"))

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness_RSA1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=median, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = medians, aes(label = size, y = fitness + 0.002)) + ylim(0.993,1.006) +
  labs(title = "TGY Hyb_t*RSA",subtitle = "Median 5roups,  down 178%")+
  theme(plot.title = element_text(size = 20))

#用errorbar-----

TGY_hybt_fitness_RSA2 <- TGY_hybt_fitness_RSA1 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(fitness = mean(fitness))
sds <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, sd)

TGY_hybt_fitness_RSA2 <- cbind(TGY_hybt_fitness_RSA2,sds)
TGY_hybt_fitness_RSA2 <- TGY_hybt_fitness_RSA2[,c(1,2,4)]
names(TGY_hybt_fitness_RSA2) <- c("group","fitness","sd")

TGY_hybt_fitness_RSA2 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean(fitness)- sd, ymax = mean(fitness)+ sd), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.01))

(means[1,2]-means[5,2])/means[1,2]

#四组-----

for (i in 1:nrow(TGY_hybt_fitness_RSA1)){
  if(TGY_hybt_fitness_RSA1[i,3]<180){
    TGY_hybt_fitness_RSA1[i,4] <-"0-180"
  }else if(TGY_hybt_fitness_RSA1[i,3]<360){
    TGY_hybt_fitness_RSA1[i,4] <-"180-360"
  }else if(TGY_hybt_fitness_RSA1[i,3]<540){
    TGY_hybt_fitness_RSA1[i,4] <-"360-540"
  }else{
    TGY_hybt_fitness_RSA1[i,4] <-"540-721"
  }}

names(TGY_hybt_fitness_RSA1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, median)%>%
  mutate(size = round(1-fitness,6))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, mean)%>%
  mutate(size = round(1-fitness,6))

TGY_hybt_fitness_RSA1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=median, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = medians, aes(label = size, y = fitness + 0.002)) + ylim(0.993,1.006) +
  labs(title = "TGY Hyb_t*RSA",subtitle = "Median 4roups,  down 151%")

(means[1,2]-means[4,2])/means[1,2]
(medians[1,2]-medians[4,2])/medians[1,2]

#(插入)判断正态性-----
shapiro.test(TGY_hybt_fitness_RSA1$fitness)
ks.test(TGY_hybt_fitness_RSA1$fitness, "pnorm", mean = mean(TGY_hybt_fitness_RSA1$fitness), sd =  sqrt(var(TGY_hybt_fitness_RSA1$fitness)))

#三组----

for (i in 1:nrow(TGY_hybt_fitness_RSA1)){
  if(TGY_hybt_fitness_RSA1[i,3]<240){
    TGY_hybt_fitness_RSA1[i,4] <-"0-240"
  }else if(TGY_hybt_fitness_RSA1[i,3]<480){
    TGY_hybt_fitness_RSA1[i,4] <-"240-480"
  }else{
    TGY_hybt_fitness_RSA1[i,4] <-"480-721"
  }}

names(TGY_hybt_fitness_RSA1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness_RSA1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = size, y = fitness + 0.002)) +
  ylim(0.993,1.006) +
  labs(title = "TGY Hyb_t*RSA",subtitle = "Mean 3roups,  down 479%")+
  theme(plot.title = element_text(size = 20))

(means[1,2]-means[3,2])/means[1,2]
(medians[1,2]-medians[3,2])/medians[1,2]

#用errorbar----

TGY_hybt_fitness_RSA2 <- TGY_hybt_fitness_RSA1 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(fitness = mean(fitness))
sds <- aggregate(fitness ~  group, TGY_hybt_fitness_RSA1, sd)

TGY_hybt_fitness_RSA2 <- cbind(TGY_hybt_fitness_RSA2,sds)
TGY_hybt_fitness_RSA2 <- TGY_hybt_fitness_RSA2[,c(1,2,4)]
names(TGY_hybt_fitness_RSA2) <- c("group","fitness","sd")

TGY_hybt_fitness_RSA2 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean(fitness)- sd, ymax = mean(fitness)+ sd), 
                width = 0.2, position = position_dodge(0.9)) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.01))+
  ylim(1,1.007)


(medians[1,2]-medians[3,2])/medians[1,2]
(means[1,2]-means[3,2])/means[1,2]


#200622 内部氨基酸-----

TGY_hybt_fitness_buried <- TGY_RSA0.2_buried[,c(11,16)]
TGY_hybt_fitness_buried1 <- TGY_hybt_fitness_buried[sort(TGY_hybt_fitness_buried$hyb_t, index.return = T)$ix,]
TGY_hybt_fitness_buried1$seq <- 1:671

TGY_hybt_fitness_buried1[,4] <- "group"

#十组-----

for (i in 1:nrow(TGY_hybt_fitness_buried1)){
  if(TGY_hybt_fitness_buried1[i,3]<66){
    TGY_hybt_fitness_buried1[i,4] <-"0-66"
  }else if(TGY_hybt_fitness_buried1[i,3]<132){
    TGY_hybt_fitness_buried1[i,4] <-"66-132"
  }else if(TGY_hybt_fitness_buried1[i,3]<198){
    TGY_hybt_fitness_buried1[i,4] <-"132-198"
  }else if(TGY_hybt_fitness_buried1[i,3]<264){
    TGY_hybt_fitness_buried1[i,4] <-"198-264"
  }else if(TGY_hybt_fitness_buried1[i,3]<330){
    TGY_hybt_fitness_buried1[i,4] <-"264-330"
  }else if(TGY_hybt_fitness_buried1[i,3]<396){
    TGY_hybt_fitness_buried1[i,4] <-"330-396"
  }else if(TGY_hybt_fitness_buried1[i,3]<462){
    TGY_hybt_fitness_buried1[i,4] <-"396-462"
  }else if(TGY_hybt_fitness_buried1[i,3]<528){
    TGY_hybt_fitness_buried1[i,4] <-"462-528"
  }else if(TGY_hybt_fitness_buried1[i,3]<594){
    TGY_hybt_fitness_buried1[i,4] <-"528-594"
  }else{
    TGY_hybt_fitness_buried1[i,4] <-"594-671"
  }}

names(TGY_hybt_fitness_buried1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness_buried1$group<- factor(TGY_hybt_fitness_buried1$group,levels=c("0-66","66-132","132-198","198-264","264-330","330-396","396-462","462-528","528-594","594-667"))

TGY_hybt_fitness_buried1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) 


#五组------

for (i in 1:nrow(TGY_hybt_fitness_buried1)){
  if(TGY_hybt_fitness_buried1[i,3]<133){
    TGY_hybt_fitness_buried1[i,4] <-"0-133"
  }else if(TGY_hybt_fitness_buried1[i,3]<266){
    TGY_hybt_fitness_buried1[i,4] <-"133-266"
  }else if(TGY_hybt_fitness_buried1[i,3]<399){
    TGY_hybt_fitness_buried1[i,4] <-"266-399"
  }else if(TGY_hybt_fitness_buried1[i,3]<532){
    TGY_hybt_fitness_buried1[i,4] <-"399-532"
  }else{
    TGY_hybt_fitness_buried1[i,4] <-"532-671"
  }}

names(TGY_hybt_fitness_buried1) <- c("hyb_t","fitness","seq","group")
TGY_hybt_fitness_buried1$group<- factor(TGY_hybt_fitness_buried1$group,levels=c("0-144","144-288","288-432","432-576","576-721"))

means <- aggregate(fitness ~  group, TGY_hybt_fitness_buried1, mean)%>%
  mutate(size = round(1-fitness,5))

TGY_hybt_fitness_buried1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = size, y = fitness + 0.005))+
  labs(title = "TGY buried Hyb_t",subtitle = "Meann 5groups")+
  theme(plot.title = element_text(size = 20)) +ylim(0.993,1.007)

#四组-----

for (i in 1:nrow(TGY_hybt_fitness_buried1)){
  if(TGY_hybt_fitness_buried1[i,3]<167){
    TGY_hybt_fitness_buried1[i,4] <-"0-167"
  }else if(TGY_hybt_fitness_buried1[i,3]<334){
    TGY_hybt_fitness_buried1[i,4] <-"167-334"
  }else if(TGY_hybt_fitness_buried1[i,3]<501){
    TGY_hybt_fitness_buried1[i,4] <-"334-501"
  }else{
    TGY_hybt_fitness_buried1[i,4] <-"501-667"
  }}

names(TGY_hybt_fitness_buried1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_buried1, median)%>%
  mutate(size = round(1-fitness,5))
means <- aggregate(fitness ~  group, TGY_hybt_fitness_buried1, mean)

TGY_hybt_fitness_buried1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = medians, aes(label = size, y = fitness + 0.002)) + ylim(0.999,1.008) +
  labs(title = "TGY buried Hyb_t",subtitle = "Median 4groups")+
  theme(plot.title = element_text(size = 20)) +ylim(0.993,1.007)

(medians[1,2]-medians[4,2])/medians[1,2]

#(插入)判断正态性-----
shapiro.test(TGY_hybt_fitness_buried1$fitness)
ks.test(TGY_hybt_fitness_buried1$fitness, "pnorm", mean = mean(TGY_hybt_fitness_buried1$fitness), sd =  sqrt(var(TGY_hybt_fitness_buried1$fitness)))

#三组-----

for (i in 1:nrow(TGY_hybt_fitness_buried1)){
  if(TGY_hybt_fitness_buried1[i,3]<223){
    TGY_hybt_fitness_buried1[i,4] <-"0-223"
  }else if(TGY_hybt_fitness_buried1[i,3]<446){
    TGY_hybt_fitness_buried1[i,4] <-"223-446"
  }else{
    TGY_hybt_fitness_buried1[i,4] <-"446-667"
  }}

names(TGY_hybt_fitness_buried1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGY_hybt_fitness_buried1, median)
means <- aggregate(fitness ~  group, TGY_hybt_fitness_buried1, mean)

TGY_hybt_fitness_buried1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.002)) + ylim(0.999,1.008) +
  labs(title = "TGY buried Hyb_t",subtitle = "Mean 3groups")

(means[1,2]-means[3,2])/means[1,2]



#看外部和内部亲水到less亲水+疏水 再做一个t.test-----


TGY_hyl_c_hyb_sur <- filter(TGY_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

sd1 <- sd(TGY_hyl_c_hyb_sur$fitness)
mean1 <- mean(TGY_hyl_c_hyb_sur$fitness)

ggplot(TGY_hyl_c_hyb_sur,aes(TGY_hyl_c_hyb_sur$hyb_t,TGY_hyl_c_hyb_sur$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  labs(title = "TGY surfaced hyl to less_hyl",subtitle = "p=0.3414, R=-0.0468")

cor.test(TGY_hyl_c_hyb_sur$hyb_t,TGY_hyl_c_hyb_sur$fitness,method = "s")


TGY_hyl_c_hyb_bur <- filter(TGY_geno, RSA<0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

ggplot(TGY_hyl_c_hyb_bur,aes(TGY_hyl_c_hyb_bur$hyb_t,TGY_hyl_c_hyb_bur$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGY_hyl_c_hyb_bur$hyb_t,TGY_hyl_c_hyb_bur$fitness,method = "s")

shapiro.test(TGY_hyl_c_hyb_bur$fitness) #两组都不是正态分布，用wilcox

wilcox.test(TGY_hyl_c_hyb_sur$fitness, TGY_hyl_c_hyb_bur$fitness)


#看外部和内部疏水到less疏水+亲水 再做一个t.test-----

TGY_hyb_c_hyl_sur <- filter(TGY_geno, RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x <hydrophobicity.y)

median2 <- median(TGY_hyb_c_hyl_sur$fitness)
my_df <- TGY_hyb_c_hyl_sur$fitness %>%
  data_frame() %>%
  mutate(dif = abs(.-median2))
mad <- median(my_df$dif)

median2 - 3*mad
median2 + 3*mad

sd2 <- sd(TGY_hyb_c_hyl_sur$fitness)
mean2 -sd2
mean2 +sd2

#& fitness<1.002 & fitness >0.997

ggplot(TGY_hyb_c_hyl_sur,aes(TGY_hyb_c_hyl_sur$hyb_t,TGY_hyb_c_hyl_sur$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  labs(title = "TGY surfaced hyb to less_hyb",subtitle = " R=0.19,P=0.092")+
  theme(plot.title = element_text(size = 20)) +ylim(0.993,1.007)

cor.test(TGY_hyb_c_hyl_sur$hyb_t,TGY_hyb_c_hyl_sur$fitness,method = "s")


TGY_hyb_c_hyl_bur <- filter(TGY_geno, RSA<0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

ggplot(TGY_hyb_c_hyl_bur,aes(TGY_hyb_c_hyl_bur$hyb_t,TGY_hyb_c_hyl_bur$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()+
  labs(title = "TGY buried hyb to less_hyb",subtitle = "p<3.3e-06, R=0.285")+
  theme(plot.title = element_text(size = 20)) 

cor.test(TGY_hyb_c_hyl_bur$hyb_t,TGY_hyb_c_hyl_bur$fitness,method = "s")

shapiro.test(TGY_hyl_c_hyb_sur$fitness) #两组都不是正态分布，用wilcox

wilcox.test(TGY_hyb_c_hyl_sur$fitness, TGY_hyb_c_hyl_bur$fitness)


# 2020-0630 fitness与wzx的进行比较-----

TGY_mut_fit <- TGY_geno[,c(5,16)]

wzx_TGY_fit <- read.csv("wzx_TGY_fitness.csv", header = F)
names(wzx_TGY_fit) <- c("mut","fitness_bef")
wzx_TGY_fit <- wzx_TGY_fit %>%
  mutate(fitness = fitness_bef ** (1/106.1))

TGY_com_fitness <- merge(TGY_mut_fit,wzx_TGY_fit, by="mut", all = T)
TGY_com_fitness <- TGY_com_fitness[,c(1,2,4)]

names(TGY_com_fitness)<- c("mut","fitness_mine","fitness_wu")
TGY_com_fitness <-na.omit(TGY_com_fitness)
  
ggplot(TGY_com_fitness,aes(TGY_com_fitness$fitness_mine,TGY_com_fitness$fitness_wu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col='red' ) + xlab("fitness_mine") + ylab("fitness_wu") +
  scale_x_continuous(breaks = seq(0.99,1.03,0.01))+
  scale_y_continuous(breaks = seq(0.99,1.03,0.01)) + style.print()

test <- filter(TGY_com_fitness, fitness_mine <1 & fitness_wu >1)

# 2020-0704 fitness与wzx的ratio进行比较-----

TGY_mut_fit <- TGY_geno[,c(5,16)]
TGY_mut_fit <- TGY_mut_fit %>%
  mutate(ratio = fitness ** (106.1))

wzx_TGY_fit <- read.csv("wzx_TGY_fitness.csv", header = F)
names(wzx_TGY_fit) <- c("mut","ratio")

TGY_com_fitness <- merge(TGY_mut_fit,wzx_TGY_fit, by="mut", all = T)
TGY_com_fitness <- TGY_com_fitness[,c(1,3,4)]

names(TGY_com_fitness)<- c("mut","ratio_mine","ratio_wu")
TGY_com_fitness <-na.omit(TGY_com_fitness)

ggplot(TGY_com_fitness,aes(TGY_com_fitness$ratio_mine,TGY_com_fitness$ratio_wu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col='red' ) + xlab("fitness_mine") + ylab("fitness_wu") +
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_y_continuous(breaks = seq(0,7,1)) + style.print()



# 2020-0630 画不同表达水平（启动子）下effect size的散点图-------------

TGY_mut_hy_fitness <- TGY_RSA0.2_surfaced[,c(5,9,10,11,16)]
TGY_mut_hy_fitness <- TGY_mut_hy_fitness %>%
  mutate(s_TGY = fitness - 1)

AGY_mut_hy_fitness <- AGY_RSA0.2_surfaced[,c(5,9,10,11,16)]
AGY_mut_hy_fitness <- AGY_mut_hy_fitness %>%
  mutate(s_AGY = fitness - 1)

TGY_AGY_hy_fit <- merge(TGY_mut_hy_fitness,AGY_mut_hy_fitness, by = c("mut","hydrophobicity.x","hydrophobicity.y","hyb_t"), all = T)
names(TGY_AGY_hy_fit) <- c("mut","hydrophobicity.x","hydrophobicity.y","hyb_t","TGY_fit","s_TGY","AGY_fit","s_AGY")
TGY_AGY_hy_fit <- na.omit(TGY_AGY_hy_fit)

TGY_AGY_hy_fit <- TGY_AGY_hy_fit %>%
  mutate(group = ifelse(hyb_t>0,"hyb","hyl"))
  
  #mutate(x_group = ifelse(s_TGY>0,"xplus0","xminux0")) %>%
  #mutate(y_group = ifelse(s_AGY>0,"yplus0","yminux0"))

#dap2b <- data.frame(ymin = -Inf, ymax = Inf,xmin = c(-0.02, 0),xmax = c(0, 0.02),period = c("S-TGY<0", "S-TGY>0"))

ggplot(TGY_AGY_hy_fit,aes(x =TGY_AGY_hy_fit$s_TGY,TGY_AGY_hy_fit$s_AGY, col = group)) +
  geom_point(size=1) + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TGY") + ylab("s_AGY") +
  style.print()

# 只要 2-4/2-5两个区域的点来画椭圆------------

test1 <- filter(TGY_AGY_hy_fit, s_AGY<0 & s_TGY<s_AGY | s_TGY>0 & s_TGY<s_AGY |s_AGY>0 & s_TGY>s_AGY)

ggplot(test1,aes(x = test1$s_TGY,test1$s_AGY)) +
  geom_point(size=1) + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TGY") + ylab("s_AGY") +
  style.print()+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = 0.018, b = 0.0025, angle = pi/7), col = '#8EE5EE', fill = '#8EE5EE', alpha = 0.02) +
  geom_ellipse(aes(x0 = -0.002, y0 = 0.002, a = 0.018, b = 0.0025, angle = pi/4),col = '#C6E2FF',fill = '#C6E2FF', alpha = 0.01) +
  coord_fixed() + annotate('text', label = '75', x = 0.017, y = 0.02, size=8, colour ='#836FFF') +
  annotate('text', label = 'Binom P=0.002', x = 0.009, y = 0.016, size=7, colour ='#000000') +
  annotate('text', label ='41', x = 0.02, y = 0.017, size=8, colour ='#7AC5CD') 

binom.test(75,116,0.5)

#d <- data.frame(x = c(-0.002,0), y = c(0.002,0), ew = c(0.02,0.02), ns = c(0.005,0.005), an = c(pi/4,pi/7))
#ggplot(d, aes(x0 = x, y0 = y, a = ew/2, b = ns/2, angle = an), fill = "blue", alpha = 0.1) +
#  geom_ellipse() +
#  coord_fixed()



ggplot(TGY_AGY_hy_fit,aes(x =TGY_AGY_hy_fit$s_TGY,TGY_AGY_hy_fit$s_AGY,col= group)) +
  geom_point(size=1) + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TGY") + ylab("s_AGY") +
  style.print() + annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "cyan") + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "blueviolet")  +
  coord_fixed() + annotate('text', label = '550', x = -0.003, y = 0.02, size=9, colour ='#000000') +
  annotate('text', label = '157', x = 0.003, y = 0.02, size=9, colour ='#000000') +
  annotate('text', label = 'Binomial< 2.2e-16', x = 0, y = 0.018, size=7, colour ='#000000') 

binom.test(550,707,0.5)

ggplot(TGY_AGY_hy_fit,aes(x =TGY_AGY_hy_fit$s_TGY,TGY_AGY_hy_fit$s_AGY,col= group)) +
  geom_point(size=1) + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TGY") + ylab("s_AGY") +
  style.print() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.1, fill = "yellow") + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1, fill = "red") +
  coord_fixed() + annotate('text', label = '264', x = 0.02, y = 0.002, size=9, colour ='#000000') +
  annotate('text', label = '443', x = 0.02, y = -0.002, size=9, colour ='#000000') +
  annotate('text', label = 'Binomial< 1.699e-11', x = 0.012, y = -0.01, size=7, colour ='#000000') 

binom.test(443,707,0.5)

display.brewer.all()

# 200710 画TGY/AGY在疏水和亲水下的fitness两组箱线图--------------

TGY_group_ef <- TGY_AGY_hy_fit[,c(6,9)]
names(TGY_group_ef) <- c('Effe_size','group')
TGY_group_ef$category <- 'TGY'

AGY_group_ef <- TGY_AGY_hy_fit[,c(8,9)]
names(AGY_group_ef) <- c('Effe_size','group')
AGY_group_ef$category <- 'AGY'

TGY_AGY_group_cate <- rbind(TGY_group_ef,AGY_group_ef)
TGY_AGY_group_cate$group<- factor(TGY_AGY_group_cate$group,levels=c("hyl","hyb"))

ggplot(TGY_AGY_group_cate, aes(category, Effe_size, fill= group)) +geom_boxplot() +
  ylim(-0.01,0.01) +
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=13), axis.title.y=element_text(size=13),legend.text=element_text(size=15)) +
  labs(title = "TGY AGY samples effect size",subtitle = expression(paste(" P =  ",1.6 %*% 10^-3)))+
  theme(plot.title = element_text(size = 20)) 

TGY_hyl <- TGY_group_ef$Effe_size[which(TGY_group_ef$group == 'hyl')]
TGY_hyb <- TGY_group_ef$Effe_size[which(TGY_group_ef$group == 'hyb')]
AGY_hyl <- AGY_group_ef$Effe_size[which(AGY_group_ef$group == 'hyl')]
AGY_hyb <- AGY_group_ef$Effe_size[which(AGY_group_ef$group == 'hyb')]

wilcox.test(TGY_hyb,TGY_hyl)
wilcox.test(AGY_hyb,AGY_hyl)
wilcox.test(TGY_group_ef$Effe_size,AGY_group_ef$Effe_size)

curve(dnorm(x,mean(TGY_hyb),sd(TGY_hyb)),xlim=c(-0.02,0.02),col="blue",lwd=3)
qqnorm(TGY_group_ef$Effe_size)
shapiro.test(TGY_group_ef$Effe_size)


# 画对角线上下------------
df_poly <- data.frame(x = c(-Inf,Inf,-Inf),y= c(-Inf, Inf, Inf))
df_poly_down <- data.frame(x = c(-Inf,Inf,Inf),y= c(-Inf, Inf, -Inf))

ggplot(TGY_AGY_hy_fit,aes(x =TGY_AGY_hy_fit$s_TGY,TGY_AGY_hy_fit$s_AGY,col=group)) +
  geom_point(size=1) + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +xlab("s_TGY") + ylab("s_AGY") +
  style.print() + geom_polygon(data = df_poly, aes(x,y), inherit.aes = FALSE, alpha = 0.2, fill = "#EE6AA7")+
  geom_polygon(data = df_poly_down, aes(x,y), inherit.aes = FALSE, alpha = 0.2, fill = "#9F79EE")+
  coord_fixed() + annotate('text', label = '400', x = 0.017, y = 0.02, size=8, colour ='#EE6AA7') +
  annotate('text', label = expression(paste(" Binom =  ",5.3 %*% 10^-4)), x = 0.009, y = 0.016, size=7, colour ='#000000') +
  annotate('text', label ='307', x = 0.02, y = 0.017, size=8, colour ='#9F79EE') 

binom.test(400,707,0.5)

# TGY_AGY 对角线上hyb和hyl----------------

test2 <- filter(TGY_AGY_hy_fit, group == 'hyb')

ggplot(test2,aes(x =test2$s_TGY,test2$s_AGY)) +
  geom_point(size=1, col = '#787878') + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +xlab("s_TGY") + ylab("s_AGY") +
  style.print() + geom_polygon(data = df_poly, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#EE6A50")+
  geom_polygon(data = df_poly_down, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#FFBC00")+
  coord_fixed() + annotate('text', label = '256', x = 0.017, y = 0.02, size=8, colour ='#EE6A50') +
  annotate('text', label = ' Binom = 0.001 ', x = 0.009, y = 0.016, size=7, colour ='#000000') +
  annotate('text', label ='186', x = 0.02, y = 0.017, size=8, colour ='#FFBC00') +
  labs(title = "TGY AGY samples effect size",subtitle = "hyb Binom =0.001")+
  theme(plot.title = element_text(size = 20)) 

binom.test(256,442,0.5)

test3 <- filter(TGY_AGY_hy_fit, group == 'hyl')

ggplot(test3,aes(x =test3$s_TGY,test3$s_AGY)) +
  geom_point(size=1, col = '#787878') + xlim(-0.02,0.02) + ylim(-0.02,0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +xlab("s_TGY") + ylab("s_AGY") +
  style.print() + geom_polygon(data = df_poly, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#EEB442")+
  geom_polygon(data = df_poly_down, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#7FFF00")+
  coord_fixed() + annotate('text', label = '144', x = 0.017, y = 0.02, size=8, colour ='#EEB442') +
  annotate('text', label = 'Binom P= 0.1764', x = 0.009, y = 0.016, size=7, colour ='#000000') +
  annotate('text', label ='121', x = 0.02, y = 0.017, size=8, colour ='#7FFF00') +
  labs(title = "TGY AGY samples effect size",subtitle = "hyl Binom =0.1764")

binom.test(144,265,0.5)


#````````````````````````````````````````
sum(TGY_AGY_hy_fit$group=='hyl')

area1 <- filter(TGY_AGY_hy_fit, s_TGY<0 & s_AGY<s_TGY)
sum(area1$group=='hyb')
sum(area1$group=='hyb')/442
sum(area1$group=='hyl')
sum(area1$group=='hyl')/265

area2 <- filter(TGY_AGY_hy_fit, s_AGY<0 & s_TGY<s_AGY)
sum(area2$group=='hyb')
sum(area2$group=='hyb')/442
sum(area2$group=='hyl')
sum(area2$group=='hyl')/265

area3 <- filter(TGY_AGY_hy_fit, s_TGY<0 & s_AGY>0)
sum(area3$group=='hyb')
sum(area3$group=='hyb')/442
sum(area3$group=='hyl')
sum(area3$group=='hyl')/265

area4 <- filter(TGY_AGY_hy_fit, s_TGY>0 & s_TGY<s_AGY)
sum(area4$group=='hyb')
sum(area4$group=='hyb')/442
sum(area4$group=='hyl')
sum(area4$group=='hyl')/265

area5 <- filter(TGY_AGY_hy_fit, s_AGY>0 & s_TGY>s_AGY)
sum(area5$group=='hyb')
sum(area5$group=='hyb')/442
sum(area5$group=='hyl')
sum(area5$group=='hyl')/265

area6 <- filter(TGY_AGY_hy_fit, s_AGY<0 & s_TGY>0)
sum(area6$group=='hyb')
sum(area6$group=='hyb')/442
sum(area6$group=='hyl')
sum(area6$group=='hyl')/265



# TGS_AGS 散点图------

TGS_mut_hy_fitness <- TGS_RSA0.2_surfaced[,c(5,9,10,11,16)]
TGS_mut_hy_fitness <- TGS_mut_hy_fitness %>%
  mutate(s_TGS = fitness - 1)

AGS_mut_hy_fitness <- AGS_RSA0.2_surfaced[,c(5,9,10,11,16)]
AGS_mut_hy_fitness <- AGS_mut_hy_fitness %>%
  mutate(s_AGS = fitness - 1)

TGS_AGS_hy_fit <- merge(TGS_mut_hy_fitness,AGS_mut_hy_fitness, by = c("mut","hydrophobicity.x","hydrophobicity.y","hyb_t"), all = T)
names(TGS_AGS_hy_fit) <- c("mut","hydrophobicity.x","hydrophobicity.y","hyb_t","TGS_fit","s_TGS","AGS_fit","s_AGS")
TGS_AGS_hy_fit <- na.omit(TGS_AGS_hy_fit)

TGS_AGS_hy_fit <- TGS_AGS_hy_fit %>%
  mutate(group = ifelse(hyb_t>0,"hyb","hyl"))

ggplot(TGS_AGS_hy_fit,aes(TGS_AGS_hy_fit$s_TGS,TGS_AGS_hy_fit$s_AGS,col = group))+
  geom_point(size=1) + xlim(-0.06,0.09) + ylim(-0.06,0.09) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("s_TGS") + ylab("s_AGS") +
  style.print()

sum(TGS_AGS_hy_fit$group=='hyb')
sum(TGS_AGS_hy_fit$group=='hyl')

area11 <- filter(TGS_AGS_hy_fit, s_TGS<0 & s_AGS<s_TGS)
sum(area11$group=='hyb')
sum(area11$group=='hyb')/442
sum(area11$group=='hyl')
sum(area11$group=='hyl')/265

area22 <- filter(TGS_AGS_hy_fit, s_AGS<0 & s_TGS<s_AGS)
sum(area22$group=='hyb')
sum(area22$group=='hyb')/442
sum(area22$group=='hyl')
sum(area22$group=='hyl')/265

area33 <- filter(TGS_AGS_hy_fit, s_TGS<0 & s_AGS>0)
sum(area33$group=='hyb')
sum(area33$group=='hyb')/442
sum(area33$group=='hyl')
sum(area33$group=='hyl')/265

area44 <- filter(TGS_AGS_hy_fit, s_TGS>0 & s_TGS<s_AGS)
sum(area44$group=='hyb')
sum(area44$group=='hyb')/442
sum(area44$group=='hyl')
sum(area44$group=='hyl')/265

area55 <- filter(TGS_AGS_hy_fit, s_AGS>0 & s_TGS>s_AGS)
sum(area55$group=='hyb')
sum(area55$group=='hyb')/442
sum(area55$group=='hyl')
sum(area55$group=='hyl')/265

area66 <- filter(TGS_AGS_hy_fit, s_AGS<0 & s_TGS>0)
sum(area66$group=='hyb')
sum(area66$group=='hyb')/442
sum(area66$group=='hyl')
sum(area66$group=='hyl')/265

# TGS_AGS 箱线图------

TGS_group_ef <- TGS_AGS_hy_fit[,c(6,9)]
names(TGS_group_ef) <- c('Effe_size','group')
TGS_group_ef$category <- 'TGS'

AGS_group_ef <- TGS_AGS_hy_fit[,c(8,9)]
names(AGS_group_ef) <- c('Effe_size','group')
AGS_group_ef$category <- 'AGS'

TGS_AGS_group_cate <- rbind(TGS_group_ef,AGS_group_ef)
TGS_AGS_group_cate$group<- factor(TGS_AGS_group_cate$group,levels=c("hyl","hyb"))

ggplot(TGS_AGS_group_cate, aes(category, Effe_size, fill= group)) +geom_boxplot() +
  ylim(-0.01,0.01) +
  theme(axis.text=element_text(size=20),axis.title.x =element_text(size=13), axis.title.y=element_text(size=13),legend.text=element_text(size=15)) +
  labs(title = "TGS AGS samples effect size",subtitle = " Wilcox P =  0.1451,      t.test P=0.0008,      Ks.test P=0.006")+
  theme(plot.title = element_text(size = 20)) 

TGS_hyl <- TGS_group_ef$Effe_size[which(TGS_group_ef$group == 'hyl')]
TGS_hyb <- TGS_group_ef$Effe_size[which(TGS_group_ef$group == 'hyb')]
AGS_hyl <- AGS_group_ef$Effe_size[which(AGS_group_ef$group == 'hyl')]
AGS_hyb <- AGS_group_ef$Effe_size[which(AGS_group_ef$group == 'hyb')]

wilcox.test(TGS_hyb,TGS_hyl)
wilcox.test(AGS_hyb,AGS_hyl)
wilcox.test(TGS_group_ef$Effe_size,AGS_group_ef$Effe_size)

t.test(TGS_hyb,TGS_hyl)
t.test(AGS_hyb,AGS_hyl)
t.test(TGS_group_ef$Effe_size,AGS_group_ef$Effe_size)

ks.test(TGS_group_ef$Effe_size,AGS_group_ef$Effe_size)
ks.test(AGS_hyb,AGS_hyl)
ks.test(TGS_hyb,TGS_hyl)

curve(dnorm(x,mean(TGS_hyb),sd(TGS_hyb)),xlim=c(-0.02,0.02),col="blue",lwd=3)
qqnorm(TGS_group_ef$Effe_size)
shapiro.test(AGS_group_ef$Effe_size)

# TGS_AGS 对角线上hyb和hyl----------------

test2 <- filter(TGS_AGS_hy_fit, group == 'hyb')

ggplot(test2,aes(x =test2$s_TGS,test2$s_AGS)) +
  geom_point(size=1, col = '#787878') + xlim(-0.07,0.07) + ylim(-0.07,0.07) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +xlab("s_TGS") + ylab("s_AGS") +
  style.print() + geom_polygon(data = df_poly, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#EE6A50")+
  geom_polygon(data = df_poly_down, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#FFBC00")+
  coord_fixed() + annotate('text', label = '230', x = 0.047, y = 0.06, size=8, colour ='#EE6A50') +
  annotate('text', label = expression(paste(" Binom =  ",1.39 %*% 10^-6)), x = 0.04, y = 0.04, size=7, colour ='#000000') +
  annotate('text', label ='137', x = 0.065, y = 0.055, size=8, colour ='#FFBC00') +
  labs(title = "TGS AGS samples effect size",subtitle = "hyb Binom =0.001")+
  theme(plot.title = element_text(size = 20)) 

binom.test(230,367,0.5)

test3 <- filter(TGS_AGS_hy_fit, group == 'hyl')

ggplot(test3,aes(x =test3$s_TGS,test3$s_AGS)) +
  geom_point(size=1, col = '#787878') + xlim(-0.07,0.07) + ylim(-0.07,0.07) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +xlab("s_TGS") + ylab("s_AGS") +
  style.print() + geom_polygon(data = df_poly, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#EEB442")+
  geom_polygon(data = df_poly_down, aes(x,y), inherit.aes = FALSE, alpha = 0.1, fill = "#7FFF00")+
  coord_fixed() + annotate('text', label = '144', x = 0.047, y = 0.06, size=8, colour ='#EEB442') +
  annotate('text', label = 'Binom P= 0.6232', x = 0.04, y = 0.04, size=7, colour ='#000000') +
  annotate('text', label ='121', x = 0.065, y = 0.055, size=8, colour ='#7FFF00') +
  labs(title = "TGS AGS samples effect size",subtitle = "hyl Binom =0.1764")

binom.test(137,265,0.5)





# 加权平均分------
p0.5<- c(96,94)
one <- c(83,86,85,87,90,84,93,90,95,95,85,85,85,98)
two <- c(90,95,94,88,73,90,85,86,95,95,98,97,96,95,95,93,99,98,95,89,86,70,95,83,96,95,89,88,90,98,85,89,90,95,97,95,96,88,95)
three <- c(79,80,64,93,91,95,85,64,66,75,95,96)
four <- c(83,90,73,95,85,77,84)
five <- c(87)
six <- c(93)

length(p0.5)*0.5 + length(one) +length(two)*2 + length(three)*3 + length(four)*4 + length(five)*5 + length(six)*6

sum(p0.5)*0.5 + sum(one) +sum(two)*2 + sum(three)*3 + sum(four)*4 + sum(five)*5 + sum(six)*6
87.85




# 200717 把disorder区域combine进来看相关性（blotplot) ------------

mut_disorder <- read.csv('mut_disorder.csv',header = F)
names(mut_disorder) <- c('mut','disorder_length')

TGY_sur_diso <- merge(TGY_geno,mut_disorder,by='mut',all = T)
TGY_sur_diso <- na.omit(TGY_sur_diso)

ggplot(TGY_sur_diso,aes(TGY_sur_diso$disorder_length,TGY_sur_diso$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso$disorder_length,TGY_sur_diso$fitness,method = "s")

TGY_sur_mean_disor <- TGY_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor,aes(TGY_sur_mean_disor$disorder_length,TGY_sur_mean_disor$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0095  rho= -0.52', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TGY IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TGY_sur_mean_disor$disorder_length,TGY_sur_mean_disor$fitness,method = "s")

IQR<-(quantile(TGY_sur_mean_disor$fitness,0.75)-quantile(TGY_sur_mean_disor$fitness,0.25))*1.5

quantile(TGY_sur_mean_disor$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor, fitness < quantile(TGY_sur_mean_disor$fitness,0.75)+IQR, fitness > quantile(TGY_sur_mean_disor$fitness,0.25)-IQR)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1459  rho= -0.32', x = 45, y = 1, size=7, colour ='#000000')
cor.test(test$disorder_length,test$fitness,method = "s")

#blotplot 10 hyb ------------

TGY_sur_diso$group <- ifelse(TGY_sur_diso$hyb_t>0,'hyb','hyl')
TGY_sur_diso_hyb <- filter(TGY_sur_diso,hyb_t >0)
TGY_sur_diso_hyl <- filter(TGY_sur_diso,hyb_t <0)

ggplot(TGY_sur_diso_hyb,aes(TGY_sur_diso_hyb$disorder_length,TGY_sur_diso_hyb$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'hyb P= 0.2622  rho= -0.052', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso_hyb$disorder_length,TGY_sur_diso_hyb$fitness,method = "s")

TGY_sur_mean_disor_hyb <- TGY_sur_diso_hyb %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_hyb,aes(TGY_sur_mean_disor_hyb$disorder_length,TGY_sur_mean_disor_hyb$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.01267  rho= -0.52', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TGY IDR hyb",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TGY_sur_mean_disor_hyb$disorder_length,TGY_sur_mean_disor_hyb$fitness,method = "s")

#不去偏离值信号更好

mean1 <-mean(TGY_sur_mean_disor_hyb$fitness)
sd1 <-sd(TGY_sur_mean_disor_hyb$fitness)

test <- filter(TGY_sur_mean_disor_hyb, fitness < mean1+1.5*sd1, fitness > mean1 -1.5*sd1)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1362  rho= -0.34', x = 45, y = 1, size=7, colour ='#000000')
cor.test(test$disorder_length,test$fitness,method = "s")

IQR<-(quantile(TGY_sur_mean_disor_hyb$fitness,0.75)-quantile(TGY_sur_mean_disor_hyb$fitness,0.25))*1.5

quantile(TGY_sur_mean_disor_hyb$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor_hyb$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor_hyb, fitness < quantile(TGY_sur_mean_disor_hyb$fitness,0.75)+IQR, fitness > quantile(TGY_sur_mean_disor_hyb$fitness,0.25)-IQR)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1362  rho= -0.34', x = 45, y = 1, size=7, colour ='#000000')
cor.test(test$disorder_length,test$fitness,method = "s")

#blotplot 10 hyl ------------

ggplot(TGY_sur_diso_hyl,aes(TGY_sur_diso_hyl$disorder_length,TGY_sur_diso_hyl$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'hyl P= 0.2935  rho= -0.0411', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso_hyl$disorder_length,TGY_sur_diso_hyl$fitness,method = "s")

TGY_sur_mean_disor_hyl <- TGY_sur_diso_hyl %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_hyl,aes(TGY_sur_mean_disor_hyl$disorder_length,TGY_sur_mean_disor_hyl$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.4  rho= -0.19', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TGY IDR hyl",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TGY_sur_mean_disor_hyl$disorder_length,TGY_sur_mean_disor_hyl$fitness,method = "s")

IQR<-(quantile(TGY_sur_mean_disor_hyl$fitness,0.75)-quantile(TGY_sur_mean_disor_hyl$fitness,0.25))*1.5

quantile(TGY_sur_mean_disor_hyl$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor_hyl$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor_hyl, fitness < quantile(TGY_sur_mean_disor_hyl$fitness,0.75)+IQR, fitness > quantile(TGY_sur_mean_disor_hyl$fitness,0.25)-IQR)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.5681  rho= -0.13', x = 45, y = 1, size=7, colour ='#000000')+
  labs(title = "TGY IDR hyl",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(test$disorder_length,test$fitness,method = "s")


# disorder wu 10-----------

mut_disorder_wu10 <- read.csv('mut_disorder_wu10.csv',header = F)
names(mut_disorder_wu10) <- c('mut','disorder_length')

TGY_sur_diso_wu10 <- merge(TGY_geno,mut_disorder_wu10,by='mut',all = T)
TGY_sur_diso_wu10 <- na.omit(TGY_sur_diso_wu10)

ggplot(TGY_sur_diso_wu10,aes(TGY_sur_diso_wu10$disorder_length,TGY_sur_diso_wu10$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.6485  rho= -0.017', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso_wu10$disorder_length,TGY_sur_diso_wu10$fitness,method = "s")

TGY_sur_mean_disor_wu10 <- TGY_sur_diso_wu10 %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_wu10,aes(TGY_sur_mean_disor_wu10$disorder_length,TGY_sur_mean_disor_wu10$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_mean_disor_wu10$disorder_length,TGY_sur_mean_disor_wu10$fitness,method = "s")
IQR<-(quantile(TGY_sur_mean_disor_wu10$fitness,0.75)-quantile(TGY_sur_mean_disor_wu10$fitness,0.25))*1.5

quantile(TGY_sur_mean_disor_wu10$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor_wu10$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor_wu10, fitness < quantile(TGY_sur_mean_disor_wu10$fitness,0.75)+IQR, fitness > quantile(TGY_sur_mean_disor_wu10$fitness,0.25)-IQR)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1459  rho= -0.32', x = 45, y = 1, size=7, colour ='#000000')
cor.test(test$disorder_length,test$fitness,method = "s")



# 200721 把disorder区域combine进来看相关性（blotplot) 20 -------------

mut_disorder_20web <- read.csv('mut_disorder_20web.csv',header = F)
names(mut_disorder_20web) <- c('mut','disorder_length')

TGY_sur_diso_20web <- merge(TGY_geno,mut_disorder_20web,by='mut',all = T)
TGY_sur_diso_20web <- na.omit(TGY_sur_diso_20web)

ggplot(TGY_sur_diso_20web,aes(TGY_sur_diso_20web$disorder_length,TGY_sur_diso_20web$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_diso_20web$disorder_length,TGY_sur_diso_20web$fitness,method = "s")

TGY_sur_mean_disor_20web <- TGY_sur_diso_20web %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_20web,aes(TGY_sur_mean_disor_20web$disorder_length,TGY_sur_mean_disor_20web$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_mean_disor_20web$disorder_length,TGY_sur_mean_disor_20web$fitness,method = "s")


IQR<-(quantile(TGY_sur_mean_disor_20web$fitness,0.75)-quantile(TGY_sur_mean_disor_20web$fitness,0.25))*1.5
quantile(TGY_sur_mean_disor_20web$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor_20web$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor_20web, fitness < quantile(TGY_sur_mean_disor_20web$fitness,0.75)+IQR, fitness > quantile(TGY_sur_mean_disor_20web$fitness,0.25)-IQR)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(test$disorder_length,test$fitness,method = "s")



# 200721 把disorder区域combine进来看相关性（blotplot) 5------------

mut_disorder_5web <- read.csv('mut_disorder_5web.csv',header = F)
names(mut_disorder_5web) <- c('mut','disorder_length')

TGY_sur_diso_5web <- merge(TGY_geno,mut_disorder_5web,by='mut',all = T)
TGY_sur_diso_5web <- na.omit(TGY_sur_diso_5web)

ggplot(TGY_sur_diso_5web,aes(TGY_sur_diso_5web$disorder_length,TGY_sur_diso_5web$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_diso_5web$disorder_length,TGY_sur_diso_5web$fitness,method = "s")

TGY_sur_mean_disor_5web <- TGY_sur_diso_5web %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_5web,aes(TGY_sur_mean_disor_5web$disorder_length,TGY_sur_mean_disor_5web$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_mean_disor_5web$disorder_length,TGY_sur_mean_disor_5web$fitness,method = "s")

# 200721 把disorder区域combine进来看相关性（blotplot) 15-------------

mut_disorder_15web <- read.csv('mut_disorder_15web.csv',header = F)
names(mut_disorder_15web) <- c('mut','disorder_length')

TGY_sur_diso_15web <- merge(TGY_geno,mut_disorder_15web,by='mut',all = T)
TGY_sur_diso_15web <- na.omit(TGY_sur_diso_15web)

ggplot(TGY_sur_diso_15web,aes(TGY_sur_diso_15web$disorder_length,TGY_sur_diso_15web$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_diso_15web$disorder_length,TGY_sur_diso_15web$fitness,method = "s")

TGY_sur_mean_disor_15web <- TGY_sur_diso_15web %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_15web,aes(TGY_sur_mean_disor_15web$disorder_length,TGY_sur_mean_disor_15web$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_mean_disor_15web$disorder_length,TGY_sur_mean_disor_15web$fitness,method = "s")


# 200720 把disorder区域combine进来看相关性（blotplot) ELM-----------------

mut_disorder_ELM <- read.csv('mut_disorder_ELM.csv',header = F)
names(mut_disorder_ELM) <- c('mut','disorder_length')

TGY_sur_diso_ELM <- merge(TGY_geno,mut_disorder_ELM,by='mut',all = T)
TGY_sur_diso_ELM <- na.omit(TGY_sur_diso_ELM)

ggplot(TGY_sur_diso_ELM,aes(TGY_sur_diso_ELM$disorder_length,TGY_sur_diso_ELM$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(TGY_sur_diso_ELM$disorder_length,TGY_sur_diso_ELM$fitness,method = "s")

TGY_sur_mean_disor_ELM <- TGY_sur_diso_ELM %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_ELM,aes(TGY_sur_mean_disor_ELM$disorder_length,TGY_sur_mean_disor_ELM$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.447  rho= -0.17', x = 50, y = 1, size=7, colour ='#000000')
cor.test(TGY_sur_mean_disor_ELM$disorder_length,TGY_sur_mean_disor_ELM$fitness,method = "s")

IQR<-(quantile(TGY_sur_mean_disor_ELM$fitness,0.75)-quantile(TGY_sur_mean_disor_ELM$fitness,0.25))*1.5
quantile(TGY_sur_mean_disor_ELM$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor_ELM$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor_ELM, fitness < quantile(TGY_sur_mean_disor_ELM$fitness,0.75)+IQR, fitness > quantile(TGY_sur_mean_disor_ELM$fitness,0.25)-IQR)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(test$disorder_length,test$fitness,method = "s")


# 200720 把disorder区域combine进来看相关性（DisEMBL) ------------------

mut_disorder_dise <- read.csv('mut_disorder_disemble.csv',header = F)
names(mut_disorder_dise) <- c('mut','disorder_length')

TGY_sur_diso_dise <- merge(TGY_geno,mut_disorder_dise,by='mut',all = T)
TGY_sur_diso_dise <- na.omit(TGY_sur_diso_dise)

ggplot(TGY_sur_diso_dise,aes(TGY_sur_diso_dise$disorder_length,TGY_sur_diso_dise$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0186  rho= 0.09', x = 60, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso_dise$disorder_length,TGY_sur_diso_dise$fitness,method = "s")

TGY_sur_mean_disor_dise <- TGY_sur_diso_dise %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_dise,aes(TGY_sur_mean_disor_dise$disorder_length,TGY_sur_mean_disor_dise$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0959  rho= 0.27', x = 30, y = 1.0005, size=7, colour ='#000000')
cor.test(TGY_sur_mean_disor_dise$disorder_length,TGY_sur_mean_disor_dise$fitness,method = "s")


# 200722 把disorder区域combine进来看相关性（IUpred2) -----

mut_disorder_iupred_short <- read.csv('mut_disorder_iupred.csv',header = F)
names(mut_disorder_iupred_short) <- c('mut','disorder_length')

TGY_sur_diso_iupred_short <- merge(TGY_geno,mut_disorder_iupred_short,by='mut',all = T)
TGY_sur_diso_iupred_short <- na.omit(TGY_sur_diso_iupred_short)

ggplot(TGY_sur_diso_iupred_short,aes(TGY_sur_diso_iupred_short$disorder_length,TGY_sur_diso_iupred_short$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1346  rho= 0.06', x = 40, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso_iupred_short$disorder_length,TGY_sur_diso_iupred_short$fitness,method = "s")

TGY_sur_mean_disor_iupred_short <- TGY_sur_diso_iupred_short %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_iupred_short,aes(TGY_sur_mean_disor_iupred_short$disorder_length,TGY_sur_mean_disor_iupred_short$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.36  rho= 0.17', x = 40, y = 1.0, size=7, colour ='#000000')
cor.test(TGY_sur_mean_disor_iupred_short$disorder_length,TGY_sur_mean_disor_iupred_short$fitness,method = "s")


# 200722 把disorder区域combine进来看相关性（IUpred2_long) ---------

mut_disorder_iupred_long <- read.csv('mut_disorder_iupred_long.csv',header = F)
names(mut_disorder_iupred_long) <- c('mut','disorder_length')

TGY_sur_diso_iupred_long <- merge(TGY_geno,mut_disorder_iupred_long,by='mut',all = T)
TGY_sur_diso_iupred_long <- na.omit(TGY_sur_diso_iupred_long)

ggplot(TGY_sur_diso_iupred_long,aes(TGY_sur_diso_iupred_long$disorder_length,TGY_sur_diso_iupred_long$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.1141  rho= 0.0589', x = 25, y = 1.01, size=7, colour ='#000000')
cor.test(TGY_sur_diso_iupred_long$disorder_length,TGY_sur_diso_iupred_long$fitness,method = "s")

TGY_sur_mean_disor_iupred_long <- TGY_sur_diso_iupred_long %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGY_sur_mean_disor_iupred_long,aes(TGY_sur_mean_disor_iupred_long$disorder_length,TGY_sur_mean_disor_iupred_long$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.57  rho= 0.12', x = 25, y = 1.001, size=7, colour ='#000000')
cor.test(TGY_sur_mean_disor_iupred_long$disorder_length,TGY_sur_mean_disor_iupred_long$fitness,method = "s")

IQR<-(quantile(TGY_sur_mean_disor_iupred_long$fitness,0.75)-quantile(TGY_sur_mean_disor_iupred_long$fitness,0.25))*1.5

quantile(TGY_sur_mean_disor_iupred_long$fitness,0.75) +IQR 
quantile(TGY_sur_mean_disor_iupred_long$fitness,0.25) -IQR

test <- filter(TGY_sur_mean_disor_iupred_long, fitness < 1.001, fitness > 0.997)
ggplot(test,aes(test$disorder_length,test$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm') + xlab("disorder_length") + ylab("fitness")+ style.print()
cor.test(test$disorder_length,test$fitness,method = "s")




