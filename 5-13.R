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

#传入竞争性培养的文件 内含有AA信息

t0_single <- read.csv("S1_single.csv",header = F)
names(t0_single) <- c("umi","DNA_pos","mut","DNA_AA","AA_pos","mut_AA")

t0_single1<- t0_single %>% 
  dplyr::group_by(umi,DNA_pos,mut,DNA_AA,AA_pos,mut_AA) %>%
  summarise(Freq0=length(mut))

#t3
t3_single <- read.csv("S3_single.csv",header = F)
names(t3_single) <- c("umi","DNA_pos","mut","DNA_AA","AA_pos","mut_AA")

t3_single1<- t3_single %>% 
  dplyr::group_by(umi,DNA_pos,mut,DNA_AA,AA_pos,mut_AA) %>%
  summarise(Freq7=length(mut))

single <- merge(t0_single1,t3_single1,by="umi",all=T)
single <- single[,-c(8:12)]
w<-which(single$mut_AA.x=="*")
single<- single[-w,]
w<-which(single$mut_AA.x=="")
single<- single[-w,]

single<-single[complete.cases(single[,7]),]  #删除第六列为NA的行
single[is.na(single)]<-0

#引入wt去离群值后的数值

s1_wt <- read.csv("S1_wt.csv",header = F)
names(s1_wt) <- c("umi_wt")
s1_wt$freq0 <- 1

s1_wt <- s1_wt %>%
  group_by(umi_wt) %>%
  dplyr::summarise(freq0=sum(freq0))

s3_wt <- read.csv("S3_wt.csv",header = F)
names(s3_wt) <- c("umi_wt")
s3_wt$freq7 <- 1

s3_wt <- s3_wt %>%
  group_by(umi_wt) %>%
  dplyr::summarise(freq7=sum(freq7))

wt_summary <- merge(s1_wt,s3_wt,by.x="umi_wt",all = T)
wt_summary[is.na(wt_summary)] <-0

wt_summary <- filter(wt_summary,freq0 >= 100)      #去除freq小于100

wt_summary <- wt_summary %>%
  mutate(wt_ratio = freq7/freq0)

a <- mean(wt_summary$wt_ratio)
b <- sd(wt_summary$wt_ratio) * 1.65

wt_summary1 <- filter(wt_summary, wt_ratio >= a-b & wt_ratio <= a+b)      #去除1.65sd离群值

wt_summary1 <- wt_summary1 %>%
  mutate(fitness = (freq7 / freq0 / 933687 * 1451962) ** (1/11.5) )

sum(wt_summary1$freq0) #933687
sum(wt_summary1$freq7) #1451962

wt_summary1$id <- row.names(wt_summary1)
ggplot(data=wt_summary1, aes(x=id,y=fitness)) +geom_point() +ylim(0,1.5)

#single去离群值 得到所需数据
single<- filter(single,Freq0 >=100)

single <- single %>%
  mutate(single_ratio = Freq7/Freq0)

a <- mean(single$single_ratio)
b <- sd(single$single_ratio) * 1.65

single1 <- filter(single, single_ratio >= a-b & single_ratio <= a+b)      #去除1.65sd离群值

single1 <- single1 %>%
  mutate(fitness = (Freq7 / Freq0 / 933687 * 1451962) ** (1/80))   #umi的fitness


#合并umi 计算每一个突变的fitness
single1 <- single1[,-9]
names(single1) <- c("umi","DNA_pos","mut","DNA_AA","AA_pos","mut_AA","Freq0","Freq7","fitness")

#single2 <- single1 %>%
#  group_by(mut,DNA_AA,AA_pos,mut_AA) %>%
#  summarise(fitness = mean(fitness))

single2 <- single1 %>%
  group_by(mut,DNA_AA,AA_pos,mut_AA) %>%
  summarise(Freq0 = sum(Freq0),Freq7 = sum(Freq7))

single2 <- single2 %>%
  mutate(fitness = (Freq7 / Freq0 / 933687 * 1451962) ** (1/80))        #每一个突变的fitness


#引入疏水性表
hyb <- read.csv("AA_meanhyb.csv",header = F, stringsAsFactors = F)
names(hyb) <- c("AA","mut_AA","hydrophobicity")
single3 <- merge(single2,hyb,by = "mut_AA",all=T)

# 氨基酸顺序和突变前AA疏水性
aa_que <- read.csv("AA_que.csv",header = F,stringsAsFactors = F,colClasses = "character")
aa_que <- t(aa_que)
aa_que <- as.data.frame(aa_que)
names(aa_que) <- c("AA_pos","AA_bef")   #没问题 aa_que是225个氨基酸顺序

aa_que$AA_pos <- as.factor(aa_que$AA_pos)%>%as.character()
aa_que$AA_bef <- as.factor(aa_que$AA_bef)%>%as.character()

hyb2 <- read.csv("AA_meanhyb.csv",header = F, stringsAsFactors = F)
names(hyb2) <- c("AA","AA_bef","hydrophobicity")

hyb2$AA_bef<-as.character(hyb2$AA_bef)

aa_que1 <- merge(aa_que,hyb2,by = "AA_bef",all=T)

aa_que1$AA_pos <- as.numeric(aa_que1$AA_pos)

# 删去顺序里没有的并merge
w<-which(single3$AA_pos=="1" | single3$AA_pos=="65" | single3$AA_pos=="66" | single3$AA_pos=="67" |single3$AA_pos>=230)
single6<- single3[-w,]

single6$AA_pos <- ifelse(single6$AA_pos <65, single6$AA_pos - 1, single6$AA_pos-4)

aa_bef_aft_hyb <- merge(single6,aa_que1,by="AA_pos",all = T)
aa_bef_aft_hyb$hydrophobicity.x <- as.numeric(aa_bef_aft_hyb$hydrophobicity.x)
aa_bef_aft_hyb$hydrophobicity.y <- as.numeric(aa_bef_aft_hyb$hydrophobicity.y)
aa_bef_aft_hyb$hyb_t <- aa_bef_aft_hyb$hydrophobicity.x - aa_bef_aft_hyb$hydrophobicity.y

#ggplot(aa_bef_aft_hyb,aes(aa_bef_aft_hyb$hyb_t,aa_bef_aft_hyb$fitness)) +
#  geom_point() +
#  geom_smooth()


# 引入ACC/RSA

ACC <- read.csv("ACC.csv",header = T, stringsAsFactors = F)
names(ACC) <- c("AA_bef", "ACC")
ACC$AA_pos <- 1:225

RSA2 <- RSA
RSA2 <- RSA2[,-1]
names(RSA2) <- c("RSA","AA_pos")

aa_bef_aft_hyb_RSA <- merge(aa_bef_aft_hyb, RSA2, by= "AA_pos", all=T)
aa_bef_aft_hyb_RSA <- merge(aa_bef_aft_hyb_RSA,ACC,by = "AA_pos", all = T)

aa_rsa_type <- aa_bef_aft_hyb_RSA
aa_rsa_type$rsa_type <- ifelse(aa_rsa_type$RSA>=0.2, "surfaced", "buried")

aa_rsa_type %>% 
  ggplot(aes(x = fitness,color=rsa_type))+
  geom_density() + style.print() 

#RSA0.2
aa_hyb_RSA0.2_surfaced <- filter(aa_bef_aft_hyb_RSA, RSA > 0.2)
aa_hyb_RSA0.2_surfaced <- filter(aa_hyb_RSA0.2_surfaced, hyb_t != 0)

aa_hyb_RSA0.2_buried <- filter(aa_bef_aft_hyb_RSA, RSA < 0.2)
aa_hyb_RSA0.2_buried <- filter(aa_hyb_RSA0.2_buried, hyb_t != 0)

wilcox.test(aa_hyb_RSA0.2_buried$fitness, aa_hyb_RSA0.2_surfaced$fitness)

aa_hyb_RSA0.2_surfaced$hyb_ACC <- aa_hyb_RSA0.2_surfaced$hyb_t * aa_hyb_RSA0.2_surfaced$ACC
ggplot(aa_hyb_RSA0.2_surfaced,aes(aa_hyb_RSA0.2_surfaced$hyb_ACC,aa_hyb_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*ACC") + ylab("fitness")+ style.print()     #横坐标为疏水性改变*ACC

cor.test(aa_hyb_RSA0.2_surfaced$hyb_ACC,aa_hyb_RSA0.2_surfaced$fitness,method = "s")

ggplot(aa_hyb_RSA0.2_surfaced,aes(aa_hyb_RSA0.2_surfaced$hyb_t,aa_hyb_RSA0.2_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(aa_hyb_RSA0.2_surfaced$hyb_t,aa_hyb_RSA0.2_surfaced$fitness,method = "s")

# 取亲水-疏水 与 疏水-没那么疏水and亲水 两组做显著差异分析

aa_RSA0.2_sur_qin_shu <- filter(aa_hyb_RSA0.2_surfaced, aa_hyb_RSA0.2_surfaced$hydrophobicity.y <0 & aa_hyb_RSA0.2_surfaced$hydrophobicity.x >0)

aa_RSA0.2_sur_shu_qin_lesshyb <- filter(aa_hyb_RSA0.2_surfaced, aa_hyb_RSA0.2_surfaced$hydrophobicity.y >0 & aa_hyb_RSA0.2_surfaced$hydrophobicity.x <aa_hyb_RSA0.2_surfaced$hydrophobicity.y)

t.test(aa_RSA0.2_sur_qin_shu$fitness, aa_RSA0.2_sur_shu_qin_lesshyb$fitness)

#取横轴大于0的
aa_hyb_RSA0.2_surfaced_x0 <- filter(aa_hyb_RSA0.2_surfaced, aa_hyb_RSA0.2_surfaced$hyb_t >0)

ggplot(aa_hyb_RSA0.2_surfaced_x0,aes(aa_hyb_RSA0.2_surfaced_x0$hyb_t,aa_hyb_RSA0.2_surfaced_x0$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(aa_hyb_RSA0.2_surfaced_x0$hyb_t,aa_hyb_RSA0.2_surfaced_x0$fitness,method = "s")


# 表面AA 只看AA改变 不看位置

aa_changed_RSA0.2_sur <- aa_hyb_RSA0.2_surfaced %>%
  group_by(AA_bef.x, mut_AA) %>%
  dplyr::summarise(hydrophobicity.x = mean(hydrophobicity.x), hydrophobicity.y = mean(hydrophobicity.y), hyb_t = mean(hyb_t)
                   , Freq7 = sum(Freq7), Freq0 = sum(Freq0), RSA = mean(RSA), ACC = mean(ACC))
aa_changed_RSA0.2_sur <- aa_changed_RSA0.2_sur %>%
  mutate(fitness = (Freq7 / Freq0 / 933687 * 1451962) ** (1/80))

ggplot(aa_changed_RSA0.2_sur,aes(aa_changed_RSA0.2_sur$hyb_t,aa_changed_RSA0.2_sur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(aa_changed_RSA0.2_sur$hyb_t,aa_changed_RSA0.2_sur$fitness,method = "s")

t.test(aa_changed_RSA0.2_sur$fitness,aa_hyb_RSA0.2_surfaced$fitness)

aa_changed_RSA0.2_sur <- aa_changed_RSA0.2_sur %>%
  mutate(hyb_t_RSA = hyb_t * RSA)

aa_changed_RSA0.2_sur <- aa_changed_RSA0.2_sur %>%
  mutate(hyb_t_ACC = hyb_t * ACC)

ggplot(aa_changed_RSA0.2_sur,aes(aa_changed_RSA0.2_sur$hyb_t_RSA,aa_changed_RSA0.2_sur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(aa_changed_RSA0.2_sur$hyb_t_RSA,aa_changed_RSA0.2_sur$fitness,method = "s")

ggplot(aa_changed_RSA0.2_sur,aes(aa_changed_RSA0.2_sur$hyb_t_ACC,aa_changed_RSA0.2_sur$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")+ style.print()

cor.test(aa_changed_RSA0.2_sur$hyb_t_ACC,aa_changed_RSA0.2_sur$fitness,method = "s")


# 横坐标改 RSA作为系数相乘
aa_hyb_RSA0.2_surfaced_1 <- aa_hyb_RSA0.2_surfaced
aa_hyb_RSA0.2_surfaced_1$hyb_t_RSA <- aa_hyb_RSA0.2_surfaced_1$hyb_t * aa_hyb_RSA0.2_surfaced_1$RSA


ggplot(aa_hyb_RSA0.2_surfaced_1,aes(aa_hyb_RSA0.2_surfaced_1$hyb_t_RSA,aa_hyb_RSA0.2_surfaced_1$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")+ style.print()

cor.test(aa_hyb_RSA0.2_surfaced_1$hyb_t_RSA,aa_hyb_RSA0.2_surfaced_1$fitness,method = "s")

#截掉了大于3小于-3的极端值
aa_hyb_RSA0.2_surfaced_2 <- filter(aa_hyb_RSA0.2_surfaced_1, hyb_t_RSA > -3 & hyb_t_RSA <3)

ggplot(aa_hyb_RSA0.2_surfaced_2,aes(aa_hyb_RSA0.2_surfaced_2$hyb_t_RSA,aa_hyb_RSA0.2_surfaced_2$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")

cor.test(aa_hyb_RSA0.2_surfaced_2$hyb_t_RSA,aa_hyb_RSA0.2_surfaced_2$fitness,method = "s")

#0.25
aa_hyb_RSA0.25_surfaced <- filter(aa_bef_aft_hyb_RSA, RSA > 0.25)
aa_hyb_RSA0.25_surfaced <- filter(aa_hyb_RSA0.25_surfaced, hyb_t != 0)

ggplot(aa_hyb_RSA0.25_surfaced,aes(aa_hyb_RSA0.25_surfaced$hyb_t,aa_hyb_RSA0.25_surfaced$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")

cor.test(aa_hyb_RSA0.25_surfaced$hyb_t,aa_hyb_RSA0.25_surfaced$fitness,method = "s")


# 0.2取绝对值

hyb_RSA0.2 <- filter(aa_hyb_RSA0.2_surfaced, aa_hyb_RSA0.2_surfaced$hyb_t >0)
hyl_RSA0.2 <- filter(aa_hyb_RSA0.2_surfaced, aa_hyb_RSA0.2_surfaced$hyb_t <0)

hyl_RSA0.2_abs <- hyl_RSA0.2
hyl_RSA0.2_abs$hyb_t <- hyl_RSA0.2_abs$hyb_t * -1

hyl_RSA0.2_abs$hyb_ACC <- hyl_RSA0.2_abs$hyb_ACC * -1

hyb_RSA0.2$type <- "hyb"
hyl_RSA0.2_abs$type <- "hyl"

hyb_hyl_abs <- rbind(hyb_RSA0.2, hyl_RSA0.2_abs)

# cutoff fig
data <- lapply(c(1:8),function(x){
  cutoff2_hyb <- filter(hyb_hyl_abs, hyb_hyl_abs$hyb_t >x & hyb_hyl_abs$type == "hyb")
cutoff2_hyl <- filter(hyb_hyl_abs, hyb_hyl_abs$hyb_t >x & hyb_hyl_abs$type == "hyl")
data1 <- data.frame(cutoff = x, hyb_mean = mean(cutoff2_hyb$fitness), hyl_mean = mean(cutoff2_hyl$fitness), hyb_med = median(cutoff2_hyb$fitness), hyl_med = median(cutoff2_hyl$fitness))
}) %>% rbind.fill()

data2 <- melt(data, id.vars = c("cutoff"))

data2_mean <- filter(data2, variable == "hyl_mean" |  variable =="hyb_mean")

ggplot(data2_mean, aes(cutoff,value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity") + coord_cartesian(ylim = c(0.99, 1)) + 
  theme(legend.position='none')+ xlab("changed hydrophobicity cutoff") + ylab("fitness mean") + style.print() 

data2_med <- filter(data2, variable == "hyl_med" |  variable =="hyb_med")

ggplot(data2_med, aes(cutoff,value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity") + coord_cartesian(ylim = c(0.99, 1)) + 
  theme(legend.position='none') + xlab("changed hydrophobicity cutoff") + ylab("fitness med") + style.print() 


ggplot(hyb_hyl_abs,aes(hyb_hyl_abs$hyb_t,hyb_hyl_abs$fitness, fill= type)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness") + style.print()   #绝对值线

wilcox.test(hyb_RSA0.2$fitness, hyl_RSA0.2_abs$fitness)

# cutoff5
hyb_RSA0.2_c5 <- filter(hyb_RSA0.2, hyb_t >=5)
hyl_RSA0.2_c5 <- filter(hyl_RSA0.2_abs, hyb_t >=5)
wilcox.test(hyb_RSA0.2_c5$fitness, hyl_RSA0.2_c5$fitness)


#RSA0.5

aa_hyb_RSA0.5_surfaced <- filter(aa_bef_aft_hyb_RSA, RSA > 0.5)
aa_hyb_RSA0.5_surfaced <- filter(aa_hyb_RSA0.5_surfaced, hyb_t != 0)

ggplot(aa_hyb_RSA0.5_surfaced,aes(aa_hyb_RSA0.5_surfaced$hyb_t,aa_hyb_RSA0.5_surfaced$fitness)) +
  geom_point() +
  geom_smooth()

cor.test(aa_hyb_RSA0.5_surfaced$hyb_t,aa_hyb_RSA0.5_surfaced$fitness,method = "s")

   #RSA0.5 乘以RSA系数

aa_hyb_RSA0.5_surfaced_1 <- aa_hyb_RSA0.5_surfaced
aa_hyb_RSA0.5_surfaced_1$hyb_t_RSA <- aa_hyb_RSA0.5_surfaced_1$hyb_t * aa_hyb_RSA0.5_surfaced_1$RSA

ggplot(aa_hyb_RSA0.5_surfaced_1,aes(aa_hyb_RSA0.5_surfaced_1$hyb_t_RSA,aa_hyb_RSA0.5_surfaced_1$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")

cor.test(aa_hyb_RSA0.5_surfaced_1$hyb_t,aa_hyb_RSA0.5_surfaced_1$fitness,method = "s")

#RSA0.36

aa_hyb_RSA0.36_surfaced <- filter(aa_bef_aft_hyb_RSA, RSA > 0.36)
aa_hyb_RSA0.36_surfaced <- filter(aa_hyb_RSA0.36_surfaced, hyb_t != 0)


ggplot(aa_hyb_RSA0.36_surfaced,aes(aa_hyb_RSA0.36_surfaced$hyb_t,aa_hyb_RSA0.36_surfaced$fitness)) +
  geom_point() +
  geom_smooth()

cor.test(aa_hyb_RSA0.36_surfaced$hyb_t,aa_hyb_RSA0.36_surfaced$fitness,method = "s")

aa_surfaced0.36_hyb_sym <- filter(aa_hyb_RSA0.36_surfaced, aa_hyb_RSA0.36_surfaced$hydrophobicity.x >0 & aa_hyb_RSA0.36_surfaced$hydrophobicity.y <0 | aa_hyb_RSA0.36_surfaced$hydrophobicity.x <0 & aa_hyb_RSA0.36_surfaced$hydrophobicity.y >0)

ggplot(aa_surfaced0.36_hyb_sym,aes(aa_surfaced0.36_hyb_sym$hyb_t,aa_surfaced0.36_hyb_sym$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")

cor.test(aa_surfaced0.36_hyb_sym$hyb_t,aa_surfaced0.36_hyb_sym$fitness,method = "s")


#取出疏水符号变化了的

aa_surfaced0.2_hyb_sym <- filter(aa_hyb_RSA0.2_surfaced, aa_hyb_RSA0.2_surfaced$hydrophobicity.x >0 & aa_hyb_RSA0.2_surfaced$hydrophobicity.y <0 | aa_hyb_RSA0.2_surfaced$hydrophobicity.x <0 & aa_hyb_RSA0.2_surfaced$hydrophobicity.y >0)

ggplot(aa_surfaced0.2_hyb_sym,aes(aa_surfaced0.2_hyb_sym$hyb_t,aa_surfaced0.2_hyb_sym$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity") + ylab("fitness")

cor.test(aa_surfaced0.2_hyb_sym$hyb_t,aa_surfaced0.2_hyb_sym$fitness,method = "s")

  #取出疏水符号变化了的 RSA0.2 乘以RSA系数

aa_surfaced0.2_hyb_sym_1 <- aa_surfaced0.2_hyb_sym
aa_surfaced0.2_hyb_sym_1$hyb_t_RSA <- aa_surfaced0.2_hyb_sym_1$hyb_t * aa_surfaced0.2_hyb_sym_1$RSA

ggplot(aa_surfaced0.2_hyb_sym_1,aes(aa_surfaced0.2_hyb_sym_1$hyb_t_RSA,aa_surfaced0.2_hyb_sym_1$fitness)) +
  geom_point() +
  geom_smooth() + xlab("changed_hydrophobicity*RSA") + ylab("fitness")

cor.test(aa_surfaced0.2_hyb_sym_1$hyb_t,aa_surfaced0.2_hyb_sym_1$fitness,method = "s")




#新疏水性表 只看最后正负疏水性
hyb_new <- read.csv("AA_hyb_new.csv",header = F)
names(hyb_new) <- c("mut_AA","hydrophobicity")
single4 <- merge(single2,hyb_new,by = "mut_AA",all=T)

single5 <- single4 %>%     #把一个密码子位点的疏水/亲水放在一起，fitness取mean
  group_by(AA_pos,hydrophobicity) %>%
  summarise(fitness = mean(fitness))

w<-which(single5$AA_pos=="1" | single5$AA_pos=="65" | single5$AA_pos=="66" | single5$AA_pos=="67" |single5$AA_pos>=230)
single5<- single5[-w,]

single5$AA_pos <- ifelse(single5$AA_pos <65, single5$AA_pos - 1, single5$AA_pos-4)

single7 <- merge(single5, RSA2, by= "AA_pos", all=T)

single7$hydrophobicity <- ifelse(single7$hydrophobicity == "a","hyb","hyl")

aa_hyb <- filter(single7,hydrophobicity == "a")
aa_hyl <- filter(single7,hydrophobicity == "b")

aa_hyb_hyl <- merge(aa_hyb,aa_hyl,by = "AA_pos",all = T)
aa_hyb_hyl <-aa_hyb_hyl[complete.cases(aa_hyb_hyl),]  #去除NA

t.test(aa_hyb_hyl$fitness.x,aa_hyb_hyl$fitness.y)

aa_hyb_hyl_surfaced0.2 <- filter(aa_hyb_hyl, RSA.x >0.2)

t.test(aa_hyb_hyl_surfaced0.2$fitness.x,aa_hyb_hyl_surfaced0.2$fitness.y)
wilcox.test(aa_hyb_hyl_surfaced0.2$fitness.x,aa_hyb_hyl_surfaced0.2$fitness.y)


single7%>%ggplot(aes(x = fitness,color=hydrophobicity))+
  geom_density()


aa_hyb_hyl <- aa_hyb_hyl %>%
  mutate(t = fitness.x - fitness.y)

pra <- aa_hyb_hyl %>%
  mutate(x = 0)




#丢弃的

#内部的
aa_hyb_RSA_buried <- filter(aa_bef_aft_hyb_RSA, RSA < 0.2)
aa_hyb_RSA_buried <- filter(aa_hyb_RSA_buried, hyb_t != 0)


ggplot(aa_hyb_RSA_buried,aes(aa_hyb_RSA_buried$hyb_t,aa_hyb_RSA_buried$fitness)) +
  geom_point() +
  geom_smooth()

cor.test(aa_hyb_RSA_buried$hyb_t,aa_hyb_RSA_buried$fitness,method = "s")

#test
ggplot(pra,aes(hyb_t,fitness)) +
  geom_point() +
  geom_smooth()

cor.test(pra$hyb_t,pra$fitness,method = "s")


# 分析蛋白内外残基

#去掉64.65.66三个位置

w<-which(single4$AA_pos=="1" | single4$AA_pos=="64" | single4$AA_pos=="65" | single4$AA_pos=="66" |single4$AA_pos>=230)
rsa_single<- single4[-w,]

rsa_single1 <- rsa_single %>%
  group_by(AA_pos) %>%
  summarise(fitness = mean(fitness))

rsa_single1 <- rsa_single1 %>%
  mutate(AA_pos = rownames(rsa_single1))

rsa_single_merge <- merge(rsa_single1,RSA1,by="AA_pos",all = T)

rsa_single_merge$AA_pos <- as.numeric(rsa_single_merge$AA_pos)

#0.2
rsa_single_merge$rsa_type <- ifelse(rsa_single_merge$RSA >= 0.2,"surfaced","buried")

rsa_single_merge %>% 
  ggplot(aes(x = fitness,color=rsa_type))+
  geom_density()

rsa_single0.2_buried <- filter(rsa_single_merge,rsa_single_merge$rsa_type == "buried")
rsa_single0.2_surfaced <- filter(rsa_single_merge,rsa_single_merge$rsa_type == "surfaced")

#0.05
rsa_single_merge$rsa_type <- ifelse(rsa_single_merge$RSA >= 0.05,"surfaced","buried")

rsa_single_merge %>% 
  ggplot(aes(x = fitness,color=rsa_type))+
  geom_density()

rsa_single0.05_buried <- filter(rsa_single_merge,rsa_single_merge$rsa_type == "buried")
rsa_single0.05_surfaced <- filter(rsa_single_merge,rsa_single_merge$rsa_type == "surfaced")

#0.3
rsa_single_merge$rsa_type <- ifelse(rsa_single_merge$RSA >= 0.3,"surfaced","buried")

rsa_single_merge %>% 
  ggplot(aes(x = fitness,color=rsa_type))+
  geom_density()

rsa_single0.3_buried <- filter(rsa_single_merge,rsa_single_merge$rsa_type == "buried")
rsa_single0.3_surfaced <- filter(rsa_single_merge,rsa_single_merge$rsa_type == "surfaced")

wilcox.test(rsa_single0.3_buried$fitness,rsa_single0.05_surfaced$fitness)
t.test(rsa_single0.3_buried$fitness,rsa_single0.05_surfaced$fitness)
