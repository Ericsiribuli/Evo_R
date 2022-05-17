AGS_day0wt <- read.csv("AGS_day0_wt.csv",header = F)
names(AGS_day0wt) <- c("umi","day0")

AGS_day7wt <- read.csv("AGS_day7_wt.csv",header = F)
names(AGS_day7wt) <- c("umi","day7")

AGS_day0single <- read.csv("AGS_day0_single.csv",header = F)
names(AGS_day0single) <- c("umi","pos","mut","DNA","AA_pos","mut_AA","day0")

AGS_day7single <- read.csv("AGS_day7_single.csv",header = F)
names(AGS_day7single) <- c("umi","pos","mut","DNA","AA_pos","mut_AA","day7")

AGS_day07_wt <- merge(AGS_day0wt,AGS_day7wt,by="umi",all = T)
AGS_day07_wt[is.na(AGS_day07_wt)]<-0
AGS_day07_wt <- AGS_day07_wt %>%
  filter(day0>0) %>%
  mutate(wt_ratio = day7/day0)

quantile(AGS_day07_wt$wt_ratio,0.975)
quantile(AGS_day07_wt$wt_ratio,0.025)
length(which((AGS_day07_wt$wt_ratio > quantile(AGS_day07_wt$wt_ratio,0.975))))
AGS_wt_ratio0 <- filter(AGS_day07_wt,wt_ratio==0)
sort(AGS_wt_ratio0$day0)[1568]
AGS_day07_wt <- filter(AGS_day07_wt,wt_ratio<quantile(AGS_day07_wt$wt_ratio,0.975) & day0 >8)

sum(AGS_day07_wt$day0)
7030535
sum(AGS_day07_wt$day7)
1469280

AGS_day07_single <- merge(AGS_day0single,AGS_day7single,by=c("umi","pos","mut","DNA","AA_pos","mut_AA"),all = T)
AGS_day07_single[is.na(AGS_day07_single)]<-0

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
AGS_day07_single <- merge(AGS_day07_single,hyb,by = "mut_AA",all=T)


# 删去顺序里没有的并merge
w<-which(AGS_day07_single$AA_pos=="1" | AGS_day07_single$AA_pos=="65" | AGS_day07_single$AA_pos=="66" | AGS_day07_single$AA_pos=="67" |AGS_day07_single$AA_pos>=230 | AGS_day07_single$mut_AA =="*")
AGS_day07_single1 <- AGS_day07_single[-w,]

AGS_day07_single1$AA_pos <- ifelse(AGS_day07_single1$AA_pos <65, AGS_day07_single1$AA_pos - 1, AGS_day07_single1$AA_pos-4)

AGS07_ba_hyb <- merge(AGS_day07_single1,aa_que1,by="AA_pos",all = T)
AGS07_ba_hyb$hydrophobicity.x <- as.numeric(AGS07_ba_hyb$hydrophobicity.x)
AGS07_ba_hyb$hydrophobicity.y <- as.numeric(AGS07_ba_hyb$hydrophobicity.y)
AGS07_ba_hyb$hyb_t <- AGS07_ba_hyb$hydrophobicity.x - AGS07_ba_hyb$hydrophobicity.y

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

AGS07_hyb_RA <- merge(AGS07_ba_hyb,ACC,by=c("AA_bef","AA_pos"),all = T)
AGS07_hyb_RA <- merge(AGS07_hyb_RA,RSA2,by=c("AA_bef","AA_pos"),all = T)

AGS_type <-AGS07_hyb_RA
AGS_type$type <- ifelse(AGS_type$RSA>=0.2, "surfaced", "buried")

AGS_type %>% 
  ggplot(aes(x = fitness,color=type))+
  geom_density() + style.print()


#200531 突变geno一样的合并

AGS_geno <- AGS07_hyb_RA %>%
  group_by(AA_bef,AA_pos,mut_AA,pos,mut,DNA,AA.x,AA.y,hydrophobicity.x,hydrophobicity.y,hyb_t,ACC,RSA) %>%
  dplyr::summarise(day0 = sum(day0),day7 = sum(day7))

AGS_geno <- filter(AGS_geno,day0 >100)

AGS_geno <- AGS_geno %>%
  mutate(fitness = (day7 / day0 / sum(AGS_day07_wt$day7) * sum(AGS_day07_wt$day0)) ** (1/74.7))


#RSA0.2
AGS_RSA0.2_surfaced <- filter(AGS_geno, RSA > 0.2 & hyb_t!=0)

AGS_RSA0.2_buried <- filter(AGS_geno, RSA < 0.2 & hyb_t != 0)

wilcox.test(AGS_RSA0.2_buried$fitness, AGS_RSA0.2_surfaced$fitness)
t.test(AGS_RSA0.2_buried$fitness, AGS_RSA0.2_surfaced$fitness)

AGS_RSA0.2_surfaced$hyb_ACC <- AGS_RSA0.2_surfaced$hyb_t * AGS_RSA0.2_surfaced$ACC
ggplot(AGS_RSA0.2_surfaced,aes(AGS_RSA0.2_surfaced$hyb_ACC,AGS_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*ACC") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*ACC

cor.test(AGS_RSA0.2_surfaced$hyb_ACC,AGS_RSA0.2_surfaced$fitness,method = "s")

ggplot(AGS_RSA0.2_surfaced,aes(AGS_RSA0.2_surfaced$hyb_t,AGS_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(AGS_RSA0.2_surfaced$hyb_t,AGS_RSA0.2_surfaced$fitness,method = "s")

AGS_RSA0.2_surfaced$hyb_RSA <- AGS_RSA0.2_surfaced$hyb_t * AGS_RSA0.2_surfaced$RSA
ggplot(AGS_RSA0.2_surfaced,aes(AGS_RSA0.2_surfaced$hyb_RSA,AGS_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*RSA

cor.test(AGS_RSA0.2_surfaced$hyb_RSA,AGS_RSA0.2_surfaced$fitness,method = "s")

# buried 0.2

ggplot(AGS_RSA0.2_buried,aes(AGS_RSA0.2_buried$hyb_t,AGS_RSA0.2_buried$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(AGS_RSA0.2_buried$hyb_t,AGS_RSA0.2_buried$fitness,method = "s")

#只看hyb>0的

AGS_0.2_surfaced_hy_more0 <- filter(AGS_geno, RSA > 0.2 & hyb_t>0)

ggplot(AGS_0.2_surfaced_hy_more0,aes(AGS_0.2_surfaced_hy_more0$hyb_t,AGS_0.2_surfaced_hy_more0$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(AGS_0.2_surfaced_hy_more0$hyb_t,AGS_0.2_surfaced_hy_more0$fitness,method = "s")

#把小于0的取绝对值，画出差距图

hyb_RSA0.2 <- filter(AGS_geno, RSA > 0.2 & hyb_t>0)
hyl_RSA0.2 <- filter(AGS_geno, RSA > 0.2 & hyb_t<0)

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

AGS_RSA0.2_cutsides4 <- filter(AGS_geno, RSA > 0.2 & hyb_t>-4 & hyb_t<4 & hyb_t!=0)

ggplot(AGS_RSA0.2_cutsides4,aes(AGS_RSA0.2_cutsides4$hyb_t,AGS_RSA0.2_cutsides4$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(AGS_RSA0.2_cutsides4$hyb_t,AGS_RSA0.2_cutsides4$fitness,method = "s")

#下列为重复（新数据除两端）
#————————————————————————————

hyb_RSA0.2 <- filter(AGS_geno, RSA > 0.2 & hyb_t<4 & hyb_t>0)
hyl_RSA0.2 <- filter(AGS_geno, RSA > 0.2 & hyb_t>-4 & hyb_t<0)

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

hyb_changed_RSA0.2 <- filter(AGS_geno, RSA>0.2 & hydrophobicity.x<0 & hydrophobicity.y>0 & hyb_t <7 & hyb_t>-7 |RSA>0.2 & hydrophobicity.x>0 & hydrophobicity.y<0 & hyb_t <7 & hyb_t>-7)

ggplot(hyb_changed_RSA0.2,aes(hyb_changed_RSA0.2$hyb_t,hyb_changed_RSA0.2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(hyb_changed_RSA0.2$hyb_t,hyb_changed_RSA0.2$fitness,method = "s")


#看亲水到疏水，疏水到没那么疏水+亲水 t.test


hyb_changedless_RSA0.2 <- filter(AGS_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y |RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

ggplot(hyb_changedless_RSA0.2,aes(hyb_changedless_RSA0.2$hyb_t,hyb_changedless_RSA0.2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(hyb_changedless_RSA0.2$hyb_t,hyb_changedless_RSA0.2$fitness,method = "s")


aa_RSA0.2_sur_qin_shu <- filter(AGS_geno, RSA>0.2 & hydrophobicity.y <0 & hydrophobicity.x > hydrophobicity.y)

aa_RSA0.2_sur_shu_qin_lesshyb <- filter(AGS_geno, RSA>0.2 & hydrophobicity.y >0 & hydrophobicity.x < hydrophobicity.y)

t.test(aa_RSA0.2_sur_qin_shu$fitness, aa_RSA0.2_sur_shu_qin_lesshyb$fitness)




#200527 把相同位置的蛋白合在一起（fitness求均值），看位置突变与fitness的关系

AGS_pos <- AGS_type %>%
  dplyr::group_by(AA_pos,AA_bef,ACC,RSA,type) %>%
  dplyr::summarise(day0 = mean(day0),day7 = mean(day7))

AGS_pos <- AGS_type %>%
  dplyr::group_by(pos,type) %>%
  dplyr::summarise(day0 = mean(day0),day7 = mean(day7))


AGS_pos <- AGS_pos %>%
  mutate(fitness = (day7 / day0 / 684106 * 1471296) ** (1/74.7))

ggplot(AGS_pos,aes(x=AA_pos,y=fitness,colour=type))+
  geom_point()


# 2020-0704 fitness与wzx的进行比较

AGS_mut_fit <- AGS_geno[,c(5,16)]

wzx_AGS_fit <- read.csv("wzx_AGS_fitness.csv", header = F)
names(wzx_AGS_fit) <- c("mut","fitness_bef")
wzx_AGS_fit <- wzx_AGS_fit %>%
  mutate(fitness = fitness_bef ** (1/74.7))

AGS_com_fitness <- merge(AGS_mut_fit,wzx_AGS_fit, by="mut", all = T)
AGS_com_fitness <- AGS_com_fitness[,c(1,2,4)]

names(AGS_com_fitness)<- c("mut","fitness_mine","fitness_wu")
AGS_com_fitness[AGS_com_fitness==0] <-NA #把有0的行换成na 下一行去掉
AGS_com_fitness <-na.omit(AGS_com_fitness)


ggplot(AGS_com_fitness,aes(AGS_com_fitness$fitness_mine,AGS_com_fitness$fitness_wu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col='red' )  + xlab("fitness_mine") + ylab("fitness_wu") +
  scale_x_continuous(breaks = seq(0.95,1.05,0.01))+
  scale_y_continuous(breaks = seq(0.95,1.05,0.01)) + style.print()

test <- filter(AGS_com_fitness, fitness_mine <1 & fitness_wu >1)


# 2020-0704 fitness与wzx的ratio进行比较

AGS_mut_fit <- AGS_geno[,c(5,16)]
AGS_mut_fit <- AGS_mut_fit %>%
  mutate(ratio = fitness ** (74.7))

wzx_AGS_fit <- read.csv("wzx_AGS_fitness.csv", header = F)
names(wzx_AGS_fit) <- c("mut","ratio")

AGS_com_fitness <- merge(AGS_mut_fit,wzx_AGS_fit, by="mut", all = T)
AGS_com_fitness <- AGS_com_fitness[,c(1,3,4)]

names(AGS_com_fitness)<- c("mut","ratio_mine","ratio_wu")
AGS_com_fitness <-na.omit(AGS_com_fitness)

AGS_com_fitness <- filter(AGS_com_fitness, ratio_mine<5&ratio_wu<5)

ggplot(AGS_com_fitness,aes(AGS_com_fitness$ratio_mine,AGS_com_fitness$ratio_wu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col='red' ) + xlab("fitness_mine") + ylab("fitness_wu") +
  scale_x_continuous(breaks = seq(0,15,1))+
  scale_y_continuous(breaks = seq(0,15,1)) + style.print()



wu_AGS_day0wt <- read.csv("wu_AGS_day0_wt.csv",header = F)
names(wu_AGS_day0wt) <- c("umi","day0")

wu_AGS_day7wt <- read.csv("wu_AGS_day7_wt.csv",header = F)
names(wu_AGS_day7wt) <- c("umi","day7")


wu_AGS_day07_wt <- merge(wu_AGS_day0wt,wu_AGS_day7wt,by="umi",all = T)
wu_AGS_day07_wt[is.na(wu_AGS_day07_wt)]<-0



