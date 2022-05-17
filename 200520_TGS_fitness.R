TGS_day0wt <- read.csv("TGS_day0_wt.csv",header = F)
names(TGS_day0wt) <- c("umi","day0")

TGS_day7wt <- read.csv("TGS_day7_wt.csv",header = F)
names(TGS_day7wt) <- c("umi","day7")

TGS_day0single <- read.csv("TGS_day0_single.csv",header = F)
names(TGS_day0single) <- c("umi","pos","mut","DNA","AA_pos","mut_AA","day0")

TGS_day7single <- read.csv("TGS_day7_single.csv",header = F)
names(TGS_day7single) <- c("umi","pos","mut","DNA","AA_pos","mut_AA","day7")

TGS_day07_wt <- merge(TGS_day0wt,TGS_day7wt,by="umi",all = T)
TGS_day07_wt[is.na(TGS_day07_wt)]<-0
TGS_day07_wt <- TGS_day07_wt %>%
  filter(day0>0) %>%
  mutate(wt_ratio = day7/day0)

quantile(TGS_day07_wt$wt_ratio,0.975)
quantile(TGS_day07_wt$wt_ratio,0.025)
length(which((TGS_day07_wt$wt_ratio > quantile(TGS_day07_wt$wt_ratio,0.975))))
TGS_wt_ratio0 <- filter(TGS_day07_wt,wt_ratio==0)
sort(TGS_wt_ratio0$day0)[1450]
TGS_day07_wt <- filter(TGS_day07_wt,wt_ratio<quantile(TGS_day07_wt$wt_ratio,0.975) & day0 >48)

sum(TGS_day07_wt$day0)
19467541
sum(TGS_day07_wt$day7)
5296667

TGS_day07_single <- merge(TGS_day0single,TGS_day7single,by=c("umi","pos","mut","DNA","AA_pos","mut_AA"), all = T)
TGS_day07_single[is.na(TGS_day07_single)]<-0

SC-74.7

# 氨基酸顺序和突变前AA疏水性
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

#引入疏水性表
hyb <- read.csv("AA_hyb.csv",header = F, stringsAsFactors = F)
names(hyb) <- c("AA","mut_AA","hydrophobicity")
TGS_day07_single <- merge(TGS_day07_single,hyb,by = "mut_AA",all=T)


# 删去顺序里没有的并merge
w<-which(TGS_day07_single$AA_pos=="1" | TGS_day07_single$AA_pos=="65" | TGS_day07_single$AA_pos=="66" | TGS_day07_single$AA_pos=="67" |TGS_day07_single$AA_pos>=230 | TGS_day07_single$mut_AA =="*")
TGS_day07_single1 <- TGS_day07_single[-w,]

TGS_day07_single1$AA_pos <- ifelse(TGS_day07_single1$AA_pos <65, TGS_day07_single1$AA_pos - 1, TGS_day07_single1$AA_pos-4)

TGS07_ba_hyb <- merge(TGS_day07_single1,aa_que1,by="AA_pos",all = T)
TGS07_ba_hyb$hydrophobicity.x <- as.numeric(TGS07_ba_hyb$hydrophobicity.x)
TGS07_ba_hyb$hydrophobicity.y <- as.numeric(TGS07_ba_hyb$hydrophobicity.y)
TGS07_ba_hyb$hyb_t <- TGS07_ba_hyb$hydrophobicity.x - TGS07_ba_hyb$hydrophobicity.y

# 引入ACC/RSA

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

# ACC/RSA merge进大数据集

TGS07_hyb_RA <- merge(TGS07_ba_hyb,ACC,by=c("AA_bef","AA_pos"),all = T)
TGS07_hyb_RA <- merge(TGS07_hyb_RA,RSA2,by=c("AA_bef","AA_pos"),all = T)

TGS_type <-TGS07_hyb_RA
TGS_type$type <- ifelse(TGS_type$RSA>=0.2, "surfaced", "buried")

TGS_type %>% 
  ggplot(aes(x = fitness,color=type))+
  geom_density() + style.print()


#200531 突变geno一样的合并

TGS_geno <- TGS07_hyb_RA %>%
  dplyr::group_by(AA_bef,AA_pos,mut_AA,pos,mut,DNA,AA.x,AA.y,hydrophobicity.x,hydrophobicity.y,hyb_t,ACC,RSA) %>%
  dplyr::summarise(day0 = sum(day0),day7 = sum(day7))

TGS_geno <- filter(TGS_geno,day0 >100)

TGS_geno <- TGS_geno %>%
  mutate(fitness = (day7 / day0 / sum(TGS_day07_wt$day7) * sum(TGS_day07_wt$day0)) ** (1/74.7))


#RSA0.2
TGS_RSA0.2_surfaced <- filter(TGS_geno, RSA > 0.2 & hyb_t!=0)

TGS_RSA0.2_buried <- filter(TGS_geno, RSA < 0.2 & hyb_t != 0)

wilcox.test(TGS_RSA0.2_buried$fitness, TGS_RSA0.2_surfaced$fitness)
t.test(TGS_RSA0.2_buried$fitness, TGS_RSA0.2_surfaced$fitness)

TGS_RSA0.2_surfaced$hyb_ACC <- TGS_RSA0.2_surfaced$hyb_t * TGS_RSA0.2_surfaced$ACC
ggplot(TGS_RSA0.2_surfaced,aes(TGS_RSA0.2_surfaced$hyb_ACC,TGS_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*ACC") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*ACC

cor.test(TGS_RSA0.2_surfaced$hyb_ACC,TGS_RSA0.2_surfaced$fitness,method = "s")

ggplot(TGS_RSA0.2_surfaced,aes(TGS_RSA0.2_surfaced$hyb_t,TGS_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_RSA0.2_surfaced$hyb_t,TGS_RSA0.2_surfaced$fitness,method = "s")

TGS_RSA0.2_surfaced$hyb_RSA <- TGS_RSA0.2_surfaced$hyb_t * TGS_RSA0.2_surfaced$RSA
ggplot(TGS_RSA0.2_surfaced,aes(TGS_RSA0.2_surfaced$hyb_RSA,TGS_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(TGS_RSA0.2_surfaced$hyb_RSA,TGS_RSA0.2_surfaced$fitness,method = "s")

# buried 0.2

ggplot(TGS_RSA0.2_buried,aes(TGS_RSA0.2_buried$hyb_t,TGS_RSA0.2_buried$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_RSA0.2_buried$hyb_t,TGS_RSA0.2_buried$fitness,method = "s")

#只看hyb>0的

TGS_0.2_surfaced_hy_more0 <- filter(TGS_geno, RSA > 0.2 & hyb_t>0)

ggplot(TGS_0.2_surfaced_hy_more0,aes(TGS_0.2_surfaced_hy_more0$hyb_t,TGS_0.2_surfaced_hy_more0$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_0.2_surfaced_hy_more0$hyb_t,TGS_0.2_surfaced_hy_more0$fitness,method = "s")

#把小于0的取绝对值，画出差距图

hyb_RSA0.2 <- filter(TGS_geno, RSA > 0.2 & hyb_t>0)
hyl_RSA0.2 <- filter(TGS_geno, RSA > 0.2 & hyb_t<0)

hyl_RSA0.2_abs <- hyl_RSA0.2
hyl_RSA0.2_abs$hyb_t <- hyl_RSA0.2_abs$hyb_t * -1

hyb_RSA0.2$type <- "hyb"
hyl_RSA0.2_abs$type <- "hyl"

hyb_hyl_abs <- rbind(hyb_RSA0.2, hyl_RSA0.2_abs)

ggplot(hyb_hyl_abs,aes(hyb_hyl_abs$hyb_t,hyb_hyl_abs$fitness, fill= type)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness") + style.print()   #绝对值线

wilcox.test(hyb_RSA0.2$fitness, hyl_RSA0.2_abs$fitness)

# 上图随着x轴增加，两条线fitness中值和平均值的差距变化

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

#把两端的值去掉（-4到4）

TGS_RSA0.2_cutsides4 <- filter(TGS_geno, RSA > 0.2 & hyb_t>-4 & hyb_t<4 & hyb_t!=0)

ggplot(TGS_RSA0.2_cutsides4,aes(TGS_RSA0.2_cutsides4$hyb_t,TGS_RSA0.2_cutsides4$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_RSA0.2_cutsides4$hyb_t,TGS_RSA0.2_cutsides4$fitness,method = "s")

#下列为重复（新数据除两端）
#————————————————————————————

hyb_RSA0.2 <- filter(TGS_geno, RSA > 0.2 & hyb_t<4 & hyb_t>0)
hyl_RSA0.2 <- filter(TGS_geno, RSA > 0.2 & hyb_t>-4 & hyb_t<0)

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


#看疏水性符号变化了的 选取-7到7

hyb_changed_RSA0.2 <- filter(TGS_geno, RSA>0.2 & hydrophobicity.x<0 & hydrophobicity.y>0 & hyb_t <7 & hyb_t>-7 |RSA>0.2 & hydrophobicity.x>0 & hydrophobicity.y<0 & hyb_t <7 & hyb_t>-7)

ggplot(hyb_changed_RSA0.2,aes(hyb_changed_RSA0.2$hyb_t,hyb_changed_RSA0.2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(hyb_changed_RSA0.2$hyb_t,hyb_changed_RSA0.2$fitness,method = "s")


#看亲水到疏水，疏水到没那么疏水+亲水 t.test


hyb_changedless_RSA0.2 <- filter(TGS_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y |RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

ggplot(hyb_changedless_RSA0.2,aes(hyb_changedless_RSA0.2$hyb_t,hyb_changedless_RSA0.2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(hyb_changedless_RSA0.2$hyb_t,hyb_changedless_RSA0.2$fitness,method = "s")


aa_RSA0.2_sur_qin_shu <- filter(TGS_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

aa_RSA0.2_sur_shu_qin_lesshyb <- filter(TGS_geno, RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

t.test(aa_RSA0.2_sur_qin_shu$fitness, aa_RSA0.2_sur_shu_qin_lesshyb$fitness)




#200527 把相同位置的蛋白合在一起（fitness求均值），看位置突变与fitness的关系

TGS_pos <- TGS_type %>%
  dplyr::group_by(AA_pos,AA_bef,ACC,RSA,type) %>%
  dplyr::summarise(day0 = mean(day0),day7 = mean(day7))

TGS_pos <- TGS_type %>%
  dplyr::group_by(pos,type) %>%
  dplyr::summarise(day0 = mean(day0),day7 = mean(day7))


TGS_pos <- TGS_pos %>%
  mutate(fitness = (day7 / day0 / sum(TGS_day07_wt$day7) * sum(TGS_day07_wt$day0)) ** (1/74.7))

ggplot(TGS_pos,aes(x=AA_pos,y=fitness,colour=type))+
  geom_point()


#200619 散点图改成分组图

TGS_hybt_fitness <- TGS_RSA0.2_surfaced[,c(11,16)]
TGS_hybt_fitness1 <- TGS_hybt_fitness[sort(TGS_hybt_fitness$hyb_t, index.return = T)$ix,]
TGS_hybt_fitness1$seq <- 1:721

TGS_hybt_fitness1[,4] <- "group"

#十组

for (i in 1:nrow(TGS_hybt_fitness1)){
  if(TGS_hybt_fitness1[i,3]<72){
    TGS_hybt_fitness1[i,4] <-"0-72"
  }else if(TGS_hybt_fitness1[i,3]<144){
    TGS_hybt_fitness1[i,4] <-"72-144"
  }else if(TGS_hybt_fitness1[i,3]<216){
    TGS_hybt_fitness1[i,4] <-"144-216"
  }else if(TGS_hybt_fitness1[i,3]<288){
    TGS_hybt_fitness1[i,4] <-"216-288"
  }else if(TGS_hybt_fitness1[i,3]<360){
    TGS_hybt_fitness1[i,4] <-"288-360"
  }else if(TGS_hybt_fitness1[i,3]<432){
    TGS_hybt_fitness1[i,4] <-"360-432"
  }else if(TGS_hybt_fitness1[i,3]<504){
    TGS_hybt_fitness1[i,4] <-"432-504"
  }else if(TGS_hybt_fitness1[i,3]<576){
    TGS_hybt_fitness1[i,4] <-"504-576"
  }else if(TGS_hybt_fitness1[i,3]<648){
    TGS_hybt_fitness1[i,4] <-"576-648"
  }else{
    TGS_hybt_fitness1[i,4] <-"648-721"
  }}

names(TGS_hybt_fitness1) <- c("hyb_t","fitness","seq","group")
TGS_hybt_fitness1$group<- factor(TGS_hybt_fitness1$group,levels=c("0-72","72-144","144-216","216-288","288-360","360-432","432-504","504-576","576-648","648-721"))

TGS_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) 


#五组

for (i in 1:nrow(TGS_hybt_fitness1)){
  if(TGS_hybt_fitness1[i,3]<144){
    TGS_hybt_fitness1[i,4] <-"0-144"
  }else if(TGS_hybt_fitness1[i,3]<288){
    TGS_hybt_fitness1[i,4] <-"144-288"
  }else if(TGS_hybt_fitness1[i,3]<432){
    TGS_hybt_fitness1[i,4] <-"288-432"
  }else if(TGS_hybt_fitness1[i,3]<576){
    TGS_hybt_fitness1[i,4] <-"432-576"
  }else{
    TGS_hybt_fitness1[i,4] <-"576-721"
  }}

names(TGS_hybt_fitness1) <- c("hyb_t","fitness","seq","group")
TGS_hybt_fitness1$group<- factor(TGS_hybt_fitness1$group,levels=c("0-144","144-288","288-432","432-576","576-721"))

means <- aggregate(fitness ~  group, TGS_hybt_fitness1, mean)

TGS_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  ylim(0.965,1.025) +
  geom_text(data = means, aes(label = fitness, y = fitness + 0.01)) +
  labs(title = "TGS Hyb_t",subtitle = "Mean 5groups,  down 0.29%")

(means[1,2]-means[5,2])/means[1,2]

#四组

for (i in 1:nrow(TGS_hybt_fitness1)){
  if(TGS_hybt_fitness1[i,3]<180){
    TGS_hybt_fitness1[i,4] <-"0-180"
  }else if(TGS_hybt_fitness1[i,3]<360){
    TGS_hybt_fitness1[i,4] <-"180-360"
  }else if(TGS_hybt_fitness1[i,3]<540){
    TGS_hybt_fitness1[i,4] <-"360-540"
  }else{
    TGS_hybt_fitness1[i,4] <-"540-721"
  }}

names(TGS_hybt_fitness1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGS_hybt_fitness1, median)
means <- aggregate(fitness ~  group, TGS_hybt_fitness1, mean)

TGS_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.02))  +
  ylim(0.965,1.025) +
  labs(title = "TGS Hyb_t",subtitle = "Mean 4groups,  down 0.16%")

(medians[1,2]-medians[4,2])/medians[1,2]
(means[1,2]-means[4,2])/means[1,2]

#(插入)判断正态性
shapiro.test(TGS_hybt_fitness1$fitness)
ks.test(TGS_hybt_fitness1$fitness, "pnorm", mean = mean(TGS_hybt_fitness1$fitness), sd =  sqrt(var(TGS_hybt_fitness1$fitness)))

#三组

for (i in 1:nrow(TGS_hybt_fitness1)){
  if(TGS_hybt_fitness1[i,3]<240){
    TGS_hybt_fitness1[i,4] <-"0-240"
  }else if(TGS_hybt_fitness1[i,3]<480){
    TGS_hybt_fitness1[i,4] <-"240-480"
  }else{
    TGS_hybt_fitness1[i,4] <-"480-721"
  }}

names(TGS_hybt_fitness1) <- c("hyb_t","fitness","seq","group")

medians <- aggregate(fitness ~  group, TGS_hybt_fitness1, median)
means <- aggregate(fitness ~  group, TGS_hybt_fitness1, mean)

TGS_hybt_fitness1 %>%
  ggplot(aes(x=group,y=fitness)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE) + 
  geom_text(data = means, aes(label = fitness, y = fitness + 0.02)) +
  ylim(0.965,1.025) +
  labs(title = "TGS Hyb_t",subtitle = "Mean 3groups,  down 0.12%")

(means[1,2]-means[3,2])/means[1,2]

#用errorbar

TGS_hybt_fitness2 <- TGS_hybt_fitness1 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(fitness = mean(fitness))
sds <- aggregate(fitness ~  group, TGS_hybt_fitness1, sd)

TGS_hybt_fitness2 <- cbind(TGS_hybt_fitness2,sds)
TGS_hybt_fitness2 <- TGS_hybt_fitness2[,c(1,2,4)]
names(TGS_hybt_fitness2) <- c("group","fitness","sd")

TGS_hybt_fitness2 %>%
  ggplot(aes(x=group,y=fitness,group=1)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean(fitness)- sd, ymax = mean(fitness)+ sd), 
                width = 0.2, position = position_dodge(0.9)) + ylim(0.95,1.05) +
  geom_text(data = means, aes(label = fitness, y = fitness + 0.04)) +
  labs(title = "TGS Hyb_t",subtitle = "Mean 3groups,  down 0.12%")




#看外部和内部亲水到疏水 再做一个t.test


TGS_hyl_c_hyb_sur <- filter(TGS_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

sd1 <- sd(TGS_hyl_c_hyb_sur$fitness)
mean1 <- mean(TGS_hyl_c_hyb_sur$fitness)

ggplot(TGS_hyl_c_hyb_sur,aes(TGS_hyl_c_hyb_sur$hyb_t,TGS_hyl_c_hyb_sur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_hyl_c_hyb_sur$hyb_t,TGS_hyl_c_hyb_sur$fitness,method = "s")


TGS_hyl_c_hyb_bur <- filter(TGS_geno, RSA<0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

ggplot(TGS_hyl_c_hyb_bur,aes(TGS_hyl_c_hyb_bur$hyb_t,TGS_hyl_c_hyb_bur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_hyl_c_hyb_bur$hyb_t,TGS_hyl_c_hyb_bur$fitness,method = "s")

shapiro.test(TGS_hyl_c_hyb_bur$fitness) #两组都不是正态分布，用wilcox

wilcox.test(TGS_hyl_c_hyb_sur$fitness, TGS_hyl_c_hyb_bur$fitness)


#看外部和内部疏水到亲水 再做一个t.test

TGS_hyb_c_hyl_sur <- filter(TGS_geno, RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x <0)

median2 <- median(TGS_hyb_c_hyl_sur$fitness)
my_df <- TGS_hyb_c_hyl_sur$fitness %>%
  data_frame() %>%
  mutate(dif = abs(.-median2))
mad <- median(my_df$dif)

median2 - 3*mad
median2 + 3*mad

sd2 <- sd(TGS_hyb_c_hyl_sur$fitness)
mean2 -sd2
mean2 +sd2

#& fitness<1.0075 & fitness >1.001

ggplot(TGS_hyb_c_hyl_sur,aes(TGS_hyb_c_hyl_sur$hyb_t,TGS_hyb_c_hyl_sur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_hyb_c_hyl_sur$hyb_t,TGS_hyb_c_hyl_sur$fitness,method = "s")


TGS_hyb_c_hyl_bur <- filter(TGS_geno, RSA<0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

ggplot(TGS_hyb_c_hyl_bur,aes(TGS_hyb_c_hyl_bur$hyb_t,TGS_hyb_c_hyl_bur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(TGS_hyb_c_hyl_bur$hyb_t,TGS_hyb_c_hyl_bur$fitness,method = "s")

shapiro.test(TGS_hyl_c_hyb_sur$fitness) #两组都不是正态分布，用wilcox

wilcox.test(TGS_hyb_c_hyl_sur$fitness, TGS_hyb_c_hyl_bur$fitness)




# 2020-0704 fitness与wzx的进行比较

TGS_mut_fit <- TGS_geno[,c(5,16)]

wzx_TGS_fit <- read.csv("wzx_TGS_fitness.csv", header = F)
names(wzx_TGS_fit) <- c("mut","fitness_bef")
wzx_TGS_fit <- wzx_TGS_fit %>%
  mutate(fitness = fitness_bef ** (1/74.7))

TGS_com_fitness <- merge(TGS_mut_fit,wzx_TGS_fit, by="mut", all = T)
TGS_com_fitness <- TGS_com_fitness[,c(1,2,4)]

names(TGS_com_fitness)<- c("mut","fitness_mine","fitness_wu")
TGS_com_fitness <-na.omit(TGS_com_fitness)

ggplot(TGS_com_fitness,aes(TGS_com_fitness$fitness_mine,TGS_com_fitness$fitness_wu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col='red' ) + xlab("fitness_mine") + ylab("fitness_wu") +
  scale_x_continuous(breaks = seq(0.95,1.05,0.01))+
  scale_y_continuous(breaks = seq(0.95,1.05,0.01)) + style.print()



#disorder----

mut_disorder <- read.csv('mut_disorder.csv',header = F)
names(mut_disorder) <- c('mut','disorder_length')

TGS_sur_diso <- merge(TGS_geno,mut_disorder,by='mut',all = T)
TGS_sur_diso <- na.omit(TGS_sur_diso)

ggplot(TGS_sur_diso,aes(TGS_sur_diso$disorder_length,TGS_sur_diso$fitness)) +
  geom_point() +
  geom_smooth() + xlab("disorder_length") + ylab("fitness")+ style.print() +
  annotate('text', label = 'P= 0.5696  rho= -0.013', x = 45, y = 1.01, size=7, colour ='#000000')
cor.test(TGS_sur_diso$disorder_length,TGS_sur_diso$fitness,method = "s")

TGS_sur_mean_disor <- TGS_sur_diso %>%
  dplyr::group_by(disorder_length) %>%
  dplyr::summarise(fitness = mean(fitness))

ggplot(TGS_sur_mean_disor,aes(TGS_sur_mean_disor$disorder_length,TGS_sur_mean_disor$fitness)) +
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("IDR length") + ylab("fitness")+ style.print()+
  annotate('text', label = 'P= 0.0095  rho= -0.52', x = 45, y = 1.0001, size=7, colour ='#000000')+
  labs(title = "TGS IDR",subtitle = "Predictor: blotplot smooth:10", size=30)
cor.test(TGS_sur_mean_disor$disorder_length,TGS_sur_mean_disor$fitness,method = "s")

