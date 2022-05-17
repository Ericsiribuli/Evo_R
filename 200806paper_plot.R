# Fig S2 ab挑选TGY  day7&day5 fitness高和低各十个做实验-----

TGY_day5_day7 <- merge(TGY_day5,TGY_day7,by='mut',all=T)
names(TGY_day5_day7) <- c('mut','fitness_day5',"fitness_day7")

ggplot(TGY_day5_day7,aes(x =TGY_day5_day7$fitness_day5,TGY_day5_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGY_fitness_day5") + ylab("TGY_fitness_day7") +
  style.print()

TGY_day5_day3 <- merge(TGY_day5,TGY_day3,by='mut',all=T)
names(TGY_day5_day3) <- c('mut','fitness_day5',"fitness_day3")

ggplot(TGY_day5_day3,aes(x =TGY_day5_day3$fitness_day3,TGY_day5_day3$fitness_day5)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGY_fitness_day3") + ylab("TGY_fitness_day5") +
  style.print()

TGY_day3_day7 <- merge(TGY_day3,TGY_day7,by='mut',all=T)
names(TGY_day3_day7) <- c('mut','fitness_day3',"fitness_day7")

ggplot(TGY_day3_day7,aes(x =TGY_day3_day7$fitness_day3,TGY_day3_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGY_fitness_day3") + ylab("TGY_fitness_day7") +
  style.print()


TGS_day5_day7 <- merge(TGS_day5,TGS_day7,by='mut',all=T)
names(TGS_day5_day7) <- c('mut','fitness_day5',"fitness_day7")
TGS_day5_day7_te <- filter(TGS_day5_day7,fitness_day5<1&fitness_day7<1)

ggplot(TGS_day5_day7,aes(x =TGS_day5_day7$fitness_day5,TGS_day5_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.95,1.1) + ylim(0.95,1.1) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGS_fitness_day5") + ylab("TGS_fitness_day7") +
  style.print()

TGY_day5_day7_te <- filter(TGY_day5_day7,fitness_day5<0.996&fitness_day7<1)

TGY_TGS_day5_day7_te <- merge(TGY_day5_day7_te,TGS_day5_day7_te,by='mut',all = F)



AGY_day5_day7 <- merge(AGY_day5,AGY_day7,by='mut',all=T)
names(AGY_day5_day7) <- c('mut','fitness_day5',"fitness_day7")

ggplot(AGY_day5_day7,aes(x =AGY_day5_day7$fitness_day5,AGY_day5_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("AGY_fitness_day5") + ylab("AGY_fitness_day7") +
  style.print()


# fitness高低的个十个点实验的生长情况-------

od_epoch2 <- read.csv('9-21epoch2_od.csv',header = T,stringsAsFactors = F)
#od_epoch2 <- as.data.frame(t(od_epoch2))
#names(od_epoch2) <- paste('t',1:120,sep = '')
od_epoch2 <- od_epoch2[,-c(1,2)]
od_epoch2_ln <-sapply(od_epoch2, function(x){log(x-0.128,exp(1))})

od_epoch_fitness <- lapply(1:96, function(x){
  doubling= 10*log(2)/ c(diff(od_epoch2_ln[,x]))
  tt <- data.frame(mean_lowfive = mean(doubling[doubling>0][order(doubling[doubling> 0])][1:5]) )
  colnames(tt)[1] <- colnames(od_epoch2_ln)[x]
  return(tt)
}) %>%cbind.data.frame()

e <- as.data.frame(t(od_epoch_fitness))
ee <- e %>% mutate(t=row.names(e)) %>% filter(e$V1 != 'NA' & e$V1 != Inf)
names(ee) <- c('doubling time','strain')
  
ee %>%
  mutate(x = factor(strain, levels = c('D2','D3','D4','E2','E3','E4','F2','F3','F4','A5','A6','A7','B5','B6','B7',
                                       'C5','C6','C7','D5','D6','D7','E5','E6','E7','F5','F6','F7','G5','G6','G7',
                                       'A8','A9','A10','B8','B9','B10',
                                       'C8','C9','C10','D8','D9','D10','E8','E9','E10','F8','F9','F10','G8','G9','G10',
                                       'B11','C11','D11','E11','F11','G11'),ordered = T)) %>%
  ggplot() +
  geom_point(aes(x=x, y= `doubling time`))


#new epoch2 data----------

od925 <- read.table('20200925TGL-SM.txt',sep = '\t',header = T,stringsAsFactors = F)
od925 <- od925[,-c(1,2)]
od925_ln <-sapply(od925,function(x){log(x-0.1251287,exp(1))}) #detract media

od925_fitness <- lapply(1:96, function(x){
  doubling= 10*log(2)/ c(diff(od925_ln[,x]))
  tt <- data.frame(mean_lowfive = mean(doubling[doubling>0][order(doubling[doubling> 0])][1:5]) )
  colnames(tt)[1] <- colnames(od925_ln)[x]
  return(tt)
}) %>%cbind.data.frame()

od925_fitness_com <- lapply(1:32, function(x){
  data <- data.frame(value=(od925_fitness[1,x*3-2]+od925_fitness[1,x*3-1]+od925_fitness[1,x*3])/3)
}) %>%rbind.fill()

od925_fitness_com$name <- c('21','2','4','6','YPD','2','4','6','0','2','4',
                            '6','1','3','5','7','1','3','5','7','1','3','5','7','1','3','5','7','2','4','6','YPD')

#928
od928 <- read.table('20200928  8-14 TGL-SM.txt',sep = '\t',header = T,stringsAsFactors = F)
od928 <- od928[,-c(1,2)]
(mean(od928$B1)+mean(od928$B2)+mean(od928$B3))/3
od928_ln <-sapply(od928,function(x){log(x-0.1262989,exp(1))}) #detract media

od928_fitness <- lapply(1:96, function(x){
  doubling= 10*log(2)/ c(diff(od928_ln[,x]))
  tt <- data.frame(mean_lowfive = mean(doubling[doubling>0][order(doubling[doubling> 0])][1:5]) )
  colnames(tt)[1] <- colnames(od928_ln)[x]
  return(tt)
}) %>%cbind.data.frame()

od928_fitness_com <- lapply(1:32, function(x){
  data <- data.frame(value=(od928_fitness[1,x*3-2]+od928_fitness[1,x*3-1]+od928_fitness[1,x*3])/3)
}) %>%rbind.fill()

od928_fitness_com$name <- c('21','9','11','13','YPD','9','11','13','0','9','11',
                            '13','8','10','12','14','8','10','12','14','8','10','12','14','8','10','12','14','9','11','13','YPD')

#1001
od1001 <- read.table('20201001.txt',sep = '\t',header = T,stringsAsFactors = F)
od1001 <- od1001[,-c(1,2)]
(mean(od1001$B1)+mean(od1001$B2)+mean(od1001$B3))/3
od1001_ln <-sapply(od1001,function(x){log(x-0.1254138,exp(1))}) #detract media

od1001_fitness <- lapply(1:96, function(x){
  doubling= 10*log(2)/ c(diff(od1001_ln[,x]))
  tt <- data.frame(mean_lowfive = mean(doubling[doubling>0][order(doubling[doubling> 0])][1:5]) )
  colnames(tt)[1] <- colnames(od1001_ln)[x]
  return(tt)
}) %>%cbind.data.frame()

od1001_fitness_com <- lapply(1:32, function(x){
  data <- data.frame(value=(od1001_fitness[1,x*3-2]+od1001_fitness[1,x*3-1]+od1001_fitness[1,x*3])/3)
}) %>%rbind.fill()

od1001_fitness_com$name <- c('21','16','18','19','YPD','16','18','19','0','16','18',
                            '19','15','16','18','na','15','17','na','19','15','17','na','na','15','17','20','20','na','17','20','20')

od_1013_all <- rbind(od925_fitness_com,od928_fitness_com,od1001_fitness_com)
od_1013_all<- na.omit(od_1013_all)
od_1013_all$name <- as.numeric(od_1013_all$name)

od_1013_all_group<- od_1013_all %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(value=mean(value))

od_1013_all_group$gene <- c('wt','A285G','G524C','A685T','G277C','A402C','G695C','C432A','T274A','T47C','T293C',
                             'T179G','T158C','C658G','A288G','G397T','C649T','G628C','G282A','G644T','C149T','BY4741')

od_1013_all_group %>%
  mutate(x = factor(gene,levels=c('wt','A285G','G524C','A685T','G277C','A402C','G695C','C432A','T274A','T47C','T293C',
                                  'T179G','T158C','C658G','A288G','G397T','C649T','G628C','G282A','G644T','C149T','BY4741'),ordered=T)) %>%
  ggplot() +
  geom_point(aes(x=x, y=value))

od_1013_all_group$fitness_od <- 57.78953/od_1013_all_group$value

od_1013_all_group$group_hl <- c('wt','high','high','high','high','high','high','high','high','high','high',
                                'low','low','low','low','low','low','low','low','low','low','by4741')

wilcox.test(od_1013_all_group$fitness_od[od_1013_all_group$group_hl=='low'],od_1013_all_group$fitness_od[od_1013_all_group$group_hl=='high'])

wilcox.test(od_1013_all_group$fitness_od[od_1013_all_group$group_hl=='low'],1)
wilcox.test(od_1013_all_group$fitness_od[od_1013_all_group$group_hl=='high'],1)


#1014 wzx改了数据

TGY_1014_fitness <- read_xlsx('10-14TGY_fitness.xlsx')
TGY_1014_fitness$group <- as.character(TGY_1014_fitness$group)

TGY_1014_fitness %>%
  ggplot() +
  geom_point(aes(x=gene,y=fitness,color=group)) +
  scale_color_manual(values = c('grey','red','blue','green'))

TGY_1014_od_fitness <- merge(TGY_1014_fitness,od_1013_all_group,all=T)

TGY_1014_od_fitness<- na.omit(TGY_1014_od_fitness)

wilcox.test(TGY_1014_od_fitness$fitness_od[TGY_1014_od_fitness$group=='low'],TGY_1014_od_fitness$fitness_od[TGY_1014_od_fitness$group=='high'])
t.test(TGY_1014_od_fitness$fitness_od[TGY_1014_od_fitness$group=='low'],TGY_1014_od_fitness$fitness_od[TGY_1014_od_fitness$group=='high'])

wilcox.test(TGY_1014_od_fitness$fitness_od[TGY_1014_od_fitness$group=='low'],1)
wilcox.test(TGY_1014_od_fitness$fitness_od[TGY_1014_od_fitness$group=='high'],1)



#Fig1 g
#TG
TG_difmut_num <- data.frame(cate=c('one mutant','two mutants','three mutants','four mutants'),num=c(2188,22127,12518,5247))

TG_difmut_num %>%
  mutate(x = factor(cate,levels=c('one mutant','two mutants','three mutants','four mutants'),ordered=T)) %>%
  ggplot(aes(x=x, y=num)) +
  geom_bar(stat="identity", fill = 'steelblue', colour = 'black', width = 0.8)+
  xlab('') + ylab('Number of captured by PacBio-CCS') +
  scale_y_continuous(expand = c(0,0),limits = c(0,25000)) +
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=15), axis.title.y=element_text(size=15)) +
  annotate('text', label = '2188', x = 1, y = 3000, size=7, colour ='#000000')+
  annotate('text', label = '22127', x = 2, y = 23000, size=7, colour ='#000000')+
  annotate('text', label = '12518', x = 3, y = 13300, size=7, colour ='#000000')+
  annotate('text', label = '5247', x = 4, y = 6000, size=7, colour ='#000000')+
  style.print()

#TU fig1g------
TU_difmut_num <- data.frame(cate=c(1,2,3,4),num=c(2375,12990,7834,3430))

p_TU_fig1g<- TU_difmut_num %>%
  mutate(x = factor(cate,levels=c(1,2,3,4),ordered=T)) %>%
  ggplot(aes(x=x, y=num)) +
  geom_bar(stat="identity", fill = 'steelblue', colour = 'black', width = 0.5)+
  xlab('Number of mutations') + ylab('Number of captured mutants') +
  scale_y_continuous(expand = c(0,0),limits = c(0,14500),breaks = seq(0,14000,2000)) +
  annotate('text', label = '2375', x = 1, y = 3300, size=6, colour ='#000000')+
  annotate('text', label = '(98.5%)', x = 1, y = 2750, size=6, colour ='#000000')+
  annotate('text', label = '12990', x = 2, y = 13320, size=6, colour ='#000000')+
  annotate('text', label = '7834', x = 3, y = 8164, size=6, colour ='#000000')+
  annotate('text', label = '3430', x = 4, y = 3770, size=6, colour ='#000000')+
  style.print() +
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18), axis.title=element_text(size=20)) 


TU_difmut_num <- data.frame(cate=c('one mutant','two mutants','three mutants','four mutants'),num=c(0.98466,0.00447121,1.681106e-06,2.450428e-09))

#p_TU_fig1g<- 
  
TU_difmut_num %>%
  mutate(x = factor(cate,levels=c('one mutant','two mutants','three mutants','four mutants'),ordered=T)) %>%
  ggplot(aes(x=x, y=num)) +
  geom_bar(stat="identity", fill = 'steelblue', colour = 'black', width = 0.5)+
  xlab('') + ylab('Fraction of captured by PacBio-CCS') +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.05)) +
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=50), axis.title.y=element_text(size=15)) +
  annotate('text', label = '98.5%', x = 1, y = 1.01, size=5, colour ='#000000')+
  annotate('text', label = '0.45%', x = 2, y = 0.03, size=5, colour ='#000000')+
  annotate('text', label = expression(paste(1.68 %*% 10^-4,' %')), x = 3, y = 0.025, size=5, colour ='#000000')+
  annotate('text', label = expression(paste(2.45 %*% 10^-7,' %')), x = 4, y = 0.02, size=5, colour ='#000000')+ 
  style.print() 


#1105 Fig1g new-----

p_TU_fig1g <- TU_difmut_num %>%
  mutate(x = factor(cate,levels=c('one mutant','two mutants','three mutants','four mutants'),ordered=T)) %>%
  ggplot(aes(x=x, y=num)) +
  #geom_bar(stat="identity", fill = 'steelblue', colour = 'black', width = 0.5)+
  xlab('') + ylab('Fraction of captured by PacBio-CCS') +
  #scale_y_continuous(expand = c(0,0),limits = c(0,1.05)) + 
  geom_segment(aes(x=x, xend=x, y=1e-9, yend=num), size=25, col='steelblue') + 
  scale_y_log10(breaks=c(0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1),expand=c(0,0),limits=c(0.000000001,10),
                labels=c(expression(paste(10^-9)),expression(paste(10^-8)),
                         expression(paste(10^-7)),expression(paste(10^-6)),expression(paste(10^-5)),
                         expression(paste(10^-4)),expression(paste(10^-3)),expression(paste(10^-2)),
                         expression(paste(10^-1)),1)) +
  annotate('text', label = '98.5%', x = 1, y = 1.5, size=5, colour ='#000000')+
  annotate('text', label = '0.45%', x = 2, y = 0.007, size=5, colour ='#000000')+
  annotate('text', label = expression(paste(1.68 %*% 10^-4,' %')), x = 3, y = 0.000003, size=5, colour ='#000000')+
  annotate('text', label = expression(paste(2.45 %*% 10^-7,' %')), x = 4, y = 0.0000000042, size=5, colour ='#000000')+ 
  style.print() + 
  theme(axis.text.x=element_text(size=19,angle = 20,vjust = 0.5),
        axis.text.y=element_text(size=17), axis.title.y=element_text(size=20))

pdf("~/fig1g_1105.pdf",width = 15,height = 14)
dev.off()

#Fig1I 相关性热图

TGY_pearson_nowt <- read.table("TGY_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(TGY_pearson_nowt) <- c("x","y","pearson")

p_TGY_nowt_pearson<-TGY_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PTDH3-GFP-YPD Repeatability of freq(No WT)",subtitle = "D1-D7 Timepoint,  S1-S3 Biological Repeat")+
  ylab('') + xlab('')


AGY_pearson_nowt <- read.table("AGY_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(AGY_pearson_nowt) <- c("x","y","pearson")

p_AGY_nowt_pearson<- AGY_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PADP1-GFP-YPD Repeatability of freq(No WT)",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+
  ylab('') + xlab('')

tiff('929scale.tiff',width = 480,height = 480)
dev.off()


#有wt的TUY
TUY_pearson <- read.csv("TUY_pearson.csv",header = F)
names(TUY_pearson) <- c("x","y","pearson")

p_TUY_pearson <- TUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradientn(colours=c("white","#FFE4E1","#FFE4E1","red"),
                       values=c(0.8,0.9899,0.99,1),
                       na.value="white", guide="colourbar",
                       name="Pearson",limits=c(0.8,1),breaks=c(0.8,0.99,1), 
                       labels=c(0.8,0.99,1))+ ylab('') + xlab('')+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size = 18))+theme(panel.border =element_blank())
  #labs(title = "TUY_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")

p2 <- p_TUY_pearson +theme(panel.border =element_blank())


pdf("~/fig1I_remove_backg.pdf",width = 15,height = 14)
dev.off()


p_TUY_pearson <- TUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low = "#FFE4E1",high = "red")+ ylab('') + xlab('')+
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size = 18))
#labs(title = "TUY_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")

pdf("~/fig1Iscale.pdf",width = 15,height = 14)
dev.off()


#10-23 new-----

p_TUY_pearson <- TUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient2(low="#FFE4E1",high="red",mid = '#FF7C74',midpoint = 0.994,limits=c(0.99,1),oob=squish,breaks=c(0.99,0.995,1), 
                       labels=c('<0.99',0.995,1)) + ylab('') + xlab('') +
  theme(axis.text.x= element_text(size=19,angle = 45,vjust = 0.9,hjust = 0.9),
        axis.text.y= element_text(size=17),text= element_text(size=20)) +theme(panel.border =element_blank())+theme(panel.grid.major=element_line(colour=NA))+
  xlab('Sample') + ylab('Sample')

pdf("~/fig1I_remove_backg1105.pdf",width = 17,height = 14)
dev.off()


pdf("~/fig1.pdf",width = 14.5,height = 7)
ggarrange(p_TU_fig1g,p_TUY_pearson,ncol=2,nrow=1,labels=c("D","E"),widths=c(0.6,0.8), heights=c(0.5,0.5),align ="hv")
dev.off()


#11-9 改成TGY的 ---------------

TG_difmut_num <- data.frame(cate=c(1,2,3,4),num=c(2183,16096,8860,3657))

p_TG_fig1g <- TG_difmut_num %>%
  mutate(x = factor(cate,levels=c(1,2,3,4),ordered=T)) %>%
  ggplot(aes(x=x, y=num)) +
  geom_bar(stat="identity", fill = 'steelblue', colour = 'black', width = 0.5)+
  xlab('Number of mutations') + ylab('Number of captured mutants') +
  scale_y_continuous(expand = c(0,0),limits = c(0,17500),breaks = seq(0,17000,2000)) +
  annotate('text', label = '2183', x = 1, y = 3450, size=6, colour ='#000000')+
  annotate('text', label = '(99.27%)', x = 1, y = 2660, size=6, colour ='#000000')+
  annotate('text', label = '16096', x = 2, y = 16596, size=6, colour ='#000000')+
  annotate('text', label = '8860', x = 3, y = 9360, size=6, colour ='#000000')+
  annotate('text', label = '3657', x = 4, y = 4150, size=6, colour ='#000000')+
  style.print() +
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18), axis.title=element_text(size=20)) 


p_TGY_pearson <- TGY_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient2(low="#FFE4E1",high="red",mid = '#FF7C74',midpoint = 0.975,limits=c(0.95,1),oob=squish,breaks=c(0.95,0.975,1), 
                       labels=c('<0.95',0.975,1)) + ylab('') + xlab('') +
  theme(axis.text.x= element_text(size=19,angle = 45,vjust = 0.9,hjust = 0.9),
        axis.text.y= element_text(size=17),text= element_text(size=20)) +theme(panel.border =element_blank())+theme(panel.grid.major=element_line(colour=NA))+
  xlab('Sample') + ylab('Sample')


pdf("~/fig1.pdf",width = 14.5,height = 7)
ggarrange(p_TG_fig1g,p_TGY_pearson,ncol=2,nrow=1,labels=c("D","E"),widths=c(0.6,0.8), heights=c(0.5,0.5),align ="hv")
dev.off()

# 1114尝试用 day0/3/7画图--------

TGY_pearson_nowt_037 <- read.table("TGY_pearson_nowt_037.txt",sep = ',', stringsAsFactors = F)
names(TGY_pearson_nowt_037) <- c("x","y","pearson")
TGY_pearson_nowt_037 <- filter(TGY_pearson_nowt_037, ! x %in% c('s1-d1','s2-d1','s3-d1'))
TGY_pearson_nowt_037 <- filter(TGY_pearson_nowt_037, ! y %in% c('s1-d1','s2-d1','s3-d1'))

#Fig1E 
P_TGY_noday5 <- TGY_pearson_nowt_037 %>%
  mutate(y = factor(y,levels=c("s1-d0","s2-d0","s1-d3","s2-d3","s3-d3","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d0","s2-d0","s1-d3","s2-d3","s3-d3","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D0-S1","D0-S2","D3-S1","D3-S2","D3-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D0-S1","D0-S2","D3-S1","D3-S2","D3-S3","D7-S1","D7-S2","Day7 - Sample3"))+
  scale_fill_gradient2(low="#FFE4E1",high="red",mid = '#FF7C74',midpoint = 0.975,limits=c(0.95,1),oob=squish,breaks=c(0.95,0.975,1), 
                       labels=c('<0.95',0.975,1)) + ylab('') + xlab('') +
  theme_bw()+ labs(fill="Pearson’s R") +
  theme(axis.text.x= element_text(size=19,angle = 45,vjust = 0.9,hjust = 0.9),
        axis.text.y= element_text(size=19),text= element_text(size=20),
        legend.position = 'top') +theme(panel.border =element_blank())+theme(panel.grid.major=element_line(colour=NA))+
  xlab('Sample') + ylab('Sample')

pdf("~/fig1.pdf",width = 13,height = 7)
ggarrange(p_TG_fig1g,P_TGY_noday5,ncol=2,nrow=1,labels=c("D","E"),widths=c(0.6,0.8), heights=c(0.5,0.5),align ="hv")
dev.off()

#AGY 热图---------

p_AGY_pearson <- AGY_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient2(low="#FFE4E1",high="red",mid = '#FF7C74',midpoint = 0.975,limits=c(0.95,1),oob=squish,breaks=c(0.95,0.975,1), 
                       labels=c('<0.95',0.975,1)) + ylab('') + xlab('') +
  theme(axis.text.x= element_text(size=19,angle = 45,vjust = 0.9,hjust = 0.9),
        axis.text.y= element_text(size=17),text= element_text(size=20)) +theme(panel.border =element_blank())+theme(panel.grid.major=element_line(colour=NA))+
  xlab('Sample') + ylab('Sample')


pdf("~/fig1.pdf",width = 14.5,height = 7)
ggarrange(p_TG_fig1g,p_AGY_pearson,ncol=2,nrow=1,labels=c("D","E"),widths=c(0.6,0.8), heights=c(0.5,0.5),align ="hv")
dev.off()

# day0/3/7画图

AGY_pearson_nowt_037 <- read.table("AGY_pearson_nowt_037.txt",sep = ',', stringsAsFactors = F)
names(AGY_pearson_nowt_037) <- c("x","y","pearson")
AGY_pearson_nowt_037 <- filter(AGY_pearson_nowt_037, ! x %in% c('s1-d1','s2-d1','s3-d1'))
AGY_pearson_nowt_037 <- filter(AGY_pearson_nowt_037, ! y %in% c('s1-d1','s2-d1','s3-d1'))

P_AGY_noday5 <- AGY_pearson_nowt_037 %>%
  mutate(y = factor(y,levels=c("s1-d0","s2-d0","s1-d3","s2-d3","s3-d3","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d0","s2-d0","s1-d3","s2-d3","s3-d3","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D0-S1","D0-S2","D3-S1","D3-S2","D3-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D0-S1","D0-S2","D3-S1","D3-S2","D3-S3","D7-S1","D7-S2","Day7 - Sample3"))+
  scale_fill_gradient2(low="#FFE4E1",high="red",mid = '#FF7C74',midpoint = 0.975,limits=c(0.95,1),oob=squish,breaks=c(0.95,0.975,1), 
                       labels=c('<0.95',0.975,1)) + ylab('') + xlab('') +
  theme_bw()+ labs(fill="Pearson’s R") +
  theme(axis.text.x= element_text(size=19,angle = 45,vjust = 0.9,hjust = 0.9),
        axis.text.y= element_text(size=19),text= element_text(size=20),
        legend.position = 'top') +theme(panel.border =element_blank())+theme(panel.grid.major=element_line(colour=NA))+
  xlab('Sample') + ylab('Sample')

pdf("~/fig1.pdf",width = 13,height = 7)
ggarrange(p_TG_fig1g,P_TGY_noday5,ncol=2,nrow=1,labels=c("D","E"),widths=c(0.6,0.8), heights=c(0.5,0.5),align ="hv")
dev.off()

-------------

#有wt的AUY
AUY_pearson <- read.csv("AUY_pearson.csv",header = F)
names(AUY_pearson) <- c("x","y","pearson")

p_AUY_pearson <-AUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AUY_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

ggarrange(p_TUY_pearson,p_AUY_pearson)

#Fig2 fitness热图-----

#TGY------
TGY_day7$mut <- as.character(TGY_day7$mut)
TGY_day7_heat <- TGY_day7

TGY_day7_heat$pos <- 
  TGY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

TGY_day7_heat$AA <- 
  TGY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

TGY_day7heat_wt <- TGY_day7_heat[,c(1,4)]
TGY_day7heat_wt$AA <- 
  TGY_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
TGY_day7heat_wt$fitness_TGY <- 1

TGY_day7_heat<- TGY_day7_heat[,c(1,2,4,5)]
TGY_day7_heat1 <- rbind(TGY_day7_heat,TGY_day7heat_wt)
TGY_day7_heat1$eff_size <- TGY_day7_heat1$fitness_TGY-1

TGY_day7_heat1$eff_size <- ifelse(TGY_day7_heat1$eff_size<=-0.01,-0.01,TGY_day7_heat1$eff_size)
TGY_day7_heat1$eff_size <- ifelse(TGY_day7_heat1$eff_size>=0.01,0.01,TGY_day7_heat1$eff_size)

Fig2_TGY <- TGY_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="#33CC00",high="#FF9900",mid='white',midpoint = 0, limits=c(-0.01,0.01),breaks=c(-0.01,0.01,0,0.005,-0.005))+
  scale_y_continuous(expand = c(0,0),limits = c(0,733))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))

#blue red


#AGY------

AGY_day7$mut <- as.character(AGY_day7$mut)
AGY_day7_heat <- AGY_day7

AGY_day7_heat$pos <- 
  AGY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

AGY_day7_heat$AA <- 
  AGY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

AGY_day7_heat <- AGY_day7_heat[,-3]
AGY_day7heat_wt <- AGY_day7_heat[,c(1,3)]
AGY_day7heat_wt$AA <- 
  AGY_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
AGY_day7heat_wt$fitness_AGY <- 1

AGY_day7_heat1 <- rbind(AGY_day7_heat,AGY_day7heat_wt)
AGY_day7_heat1$eff_size <- AGY_day7_heat1$fitness_AGY-1

AGY_day7_heat1$eff_size <- ifelse(AGY_day7_heat1$eff_size<=-0.01,-0.01,AGY_day7_heat1$eff_size)
AGY_day7_heat1$eff_size <- ifelse(AGY_day7_heat1$eff_size>=0.01,0.01,AGY_day7_heat1$eff_size)

Fig2_AGY <- AGY_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="#33CC00",high="#FF9900",mid='white',midpoint = 0, limits=c(-0.01,0.01),breaks=c(-0.01,0.01,0,0.005,-0.005))+
  scale_y_continuous(expand = c(0,0),limits = c(0,733))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))

pdf("~/fig2_GFP.pdf",width = 7,height = 20)
ggarrange(Fig2_TGY,Fig2_AGY)
dev.off()




#TUY------
TUY_day7$mut <- as.character(TUY_day7$mut)
TUY_day7_heat <- TUY_day7

TUY_day7_heat$pos <- 
  TUY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

TUY_day7_heat$AA <- 
  TUY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

TUY_day7_heat <- TUY_day7_heat[,-c(3,4)]
TUY_day7heat_wt <- TUY_day7_heat[,c(1,3)]
TUY_day7heat_wt$AA <- 
  TUY_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
TUY_day7heat_wt$fitness_TUY <- 1

TUY_day7_heat1 <- rbind(TUY_day7_heat,TUY_day7heat_wt)
TUY_day7_heat1$eff_size <- TUY_day7_heat1$fitness_TUY-1
TUY_day7_heat1$eff_size <- ifelse(TUY_day7_heat1$eff_size<=-0.03,-0.03,TUY_day7_heat1$eff_size)
TUY_day7_heat1$eff_size <- ifelse(TUY_day7_heat1$eff_size>=0.03,0.03,TUY_day7_heat1$eff_size)
  
Fig2_TUY <- TUY_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="blue",high="red",mid='white',midpoint = 0)+
  scale_y_continuous(expand = c(0,0),limits = c(0,804))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))
#blue red

#AUY-----
AUY_day7$mut <- as.character(AUY_day7$mut)
AUY_day7_heat <- AUY_day7

AUY_day7_heat$pos <- 
  AUY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

AUY_day7_heat$AA <- 
  AUY_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

AUY_day7_heat <- AUY_day7_heat[,-c(3,4)]
AUY_day7heat_wt <- AUY_day7_heat[,c(1,3)]
AUY_day7heat_wt$AA <- 
  AUY_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
AUY_day7heat_wt$fitness_AUY <- 1
AUY_day7heat_wt$pos <- 
  AUY_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

AUY_day7_heat1 <- rbind(AUY_day7_heat,AUY_day7heat_wt)
AUY_day7_heat1$eff_size <- AUY_day7_heat1$fitness_AUY-1
AUY_day7_heat1$eff_size <- ifelse(AUY_day7_heat1$eff_size<=-0.03,-0.03,AUY_day7_heat1$eff_size)
AUY_day7_heat1$eff_size <- ifelse(AUY_day7_heat1$eff_size>=0.03,0.03,AUY_day7_heat1$eff_size)

Fig2_AUY <- AUY_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="blue",high="red",mid='white',midpoint = 0,
                       limits=c(-0.03,0.03),breaks=c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03), 
                       labels=c(-0.03,-0.02,-0.01,0,0.01,0.02,0.03))+
  scale_y_continuous(expand = c(0,0),limits = c(0,804))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))+
  scale_colour_gradient(breaks = c(-0.01,0.05))

pdf("~/fig2_URA3.pdf",width = 7,height = 20)
ggarrange(Fig2_TUY,Fig2_AUY)
dev.off()

#TUU--------
TUU_day7$mut <- as.character(TUU_day7$mut)
TUU_day7_heat <- TUU_day7

TUU_day7_heat$pos <- 
  TUU_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

TUU_day7_heat$AA <- 
  TUU_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

TUU_day7_heat <- TUU_day7_heat[,-c(3)]
TUU_day7heat_wt <- TUU_day7_heat[,c(1,3)]
TUU_day7heat_wt$AA <- 
  TUU_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
TUU_day7heat_wt$fitness_TUU <- 1

TUU_day7_heat1 <- rbind(TUU_day7_heat,TUU_day7heat_wt)
TUU_day7_heat1$eff_size <- TUU_day7_heat1$fitness_TUU-1
TUU_day7_heat1$eff_size <- ifelse(TUU_day7_heat1$eff_size<=-0.06,-0.06,TUU_day7_heat1$eff_size)
TUU_day7_heat1$eff_size <- ifelse(TUU_day7_heat1$eff_size>=0.06,0.06,TUU_day7_heat1$eff_size)

Fig2_TUU <- TUU_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="blue",high="red",mid='white',midpoint = 0)+
  scale_y_continuous(expand = c(0,0),limits = c(0,804))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))
#blue red


#AUU-----
AUU_day7$mut <- as.character(AUU_day7$mut)
AUU_day7_heat <- AUU_day7

AUU_day7_heat$pos <- 
  AUU_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

AUU_day7_heat$AA <- 
  AUU_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

AUU_day7_heat <- AUU_day7_heat[,-3]
AUU_day7heat_wt <- AUU_day7_heat[,c(1,3)]
AUU_day7heat_wt$AA <- 
  AUU_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
AUU_day7heat_wt$fitness_AUU <- 1

AUU_day7_heat1 <- rbind(AUU_day7_heat,AUU_day7heat_wt)
AUU_day7_heat1$eff_size <- AUU_day7_heat1$fitness_AUU-1
AUU_day7_heat1$eff_size <- ifelse(AUU_day7_heat1$eff_size<=-0.06,-0.06,AUU_day7_heat1$eff_size)
AUU_day7_heat1$eff_size <- ifelse(AUU_day7_heat1$eff_size>=0.06,0.06,AUU_day7_heat1$eff_size)


Fig2_AUU <- AUU_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="blue",high="red",mid='white',midpoint = 0)+
  scale_y_continuous(expand = c(0,0),limits = c(0,804))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))

pdf("~/fig2_UU.pdf",width = 7,height = 20)
ggarrange(Fig2_TUU,Fig2_AUU)
dev.off()


#TUS--------
TUS_day7$mut <- as.character(TUS_day7$mut)
TUS_day7_heat <- TUS_day7

TUS_day7_heat$pos <- 
  TUS_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

TUS_day7_heat$AA <- 
  TUS_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

TUS_day7_heat <- TUS_day7_heat[,-c(3,4)]
TUS_day7heat_wt <- TUS_day7_heat[,c(1,3)]
TUS_day7heat_wt$AA <- 
  TUS_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
TUS_day7heat_wt$fitness_TUS <- 1

TUS_day7_heat1 <- rbind(TUS_day7_heat,TUS_day7heat_wt)
TUS_day7_heat1$eff_size <- TUS_day7_heat1$fitness_TUS-1
TUS_day7_heat1$eff_size <- ifelse(TUS_day7_heat1$eff_size<=-0.05,-0.05,TUS_day7_heat1$eff_size)
TUS_day7_heat1$eff_size <- ifelse(TUS_day7_heat1$eff_size>=0.05,0.05,TUS_day7_heat1$eff_size)

Fig2_TUS <- TUS_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="blue",high="red",mid='white',midpoint = 0)+
  scale_y_continuous(expand = c(0,0),limits = c(0,804))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))
#blue red


#AUS-----
AUS_day7$mut <- as.character(AUS_day7$mut)
AUS_day7_heat <- AUS_day7

AUS_day7_heat$pos <- 
  AUS_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){as.numeric(paste(x[2:4],collapse="",sep=""))}else if(length(x)==4)
  {as.numeric(paste(x[2:3],collapse="",sep=""))}else if(length(x)==3){as.numeric(paste(x[2]))}}) %>% 
  unlist()

AUS_day7_heat$AA <- 
  AUS_day7_heat$mut %>%
  strsplit("") %>% 
  lapply(function(x){if(length(x)==5){paste(x[5],collapse="",sep="")}else if(length(x)==4)
  {paste(x[4],collapse="",sep="")}else if(length(x)==3){paste(x[3],collapse="",sep="")}}) %>% 
  unlist()

AUS_day7_heat <- AUS_day7_heat[,-c(3,4)]
AUS_day7heat_wt <- AUS_day7_heat[,c(1,3)]
AUS_day7heat_wt$AA <- 
  AUS_day7heat_wt$mut %>%
  strsplit("") %>% 
  lapply(function(x){x[1]}) %>% 
  unlist()
AUS_day7heat_wt$fitness_AUS <- 1

AUS_day7_heat1 <- rbind(AUS_day7_heat,AUS_day7heat_wt)
AUS_day7_heat1$eff_size <- AUS_day7_heat1$fitness_AUS-1
AUS_day7_heat1$eff_size <- ifelse(AUS_day7_heat1$eff_size<=-0.05,-0.05,AUS_day7_heat1$eff_size)
AUS_day7_heat1$eff_size <- ifelse(AUS_day7_heat1$eff_size>=0.05,0.05,AUS_day7_heat1$eff_size)


Fig2_AUS <- AUS_day7_heat1 %>%
  mutate(mutInto = factor(AA,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutInto,y=pos,fill=eff_size )) +
  geom_tile() + xlab('') + ylab('') +
  scale_fill_gradient2(low="blue",high="red",mid='white',midpoint = 0)+
  scale_y_continuous(expand = c(0,0),limits = c(0,804))+
  theme_bw() + labs(fill="Fitness effect") +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.text=element_text(size=10),
        legend.position = 'top',
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.1))


pdf("~/fig2_US.pdf",width = 7,height = 20)
ggarrange(Fig2_TUS,Fig2_AUS)
dev.off()




pdf("~/fig2_fitness.pdf",width = 28,height = 15)
ggarrange(Fig2_TGY,Fig2_AGY,Fig2_TUY,Fig2_AUY,Fig2_TUS,Fig2_AUS,Fig2_TUU,Fig2_AUU,nrow = 1,ncol = 8,
          labels=c("A","B",'C',"D","E",'F',"G","H"))
dev.off()

pdf("~/fig2_qq.pdf",width = 150,height = 20)
ggarrange(p_TGY_sy_qq,p_AGY_sy_qq,p_TUY_sy_qq,p_AUY_sy_qq,p_TUS_sy_qq,p_AUS_sy_qq,
          p_TUU_sy_qq,p_AUU_sy_qq,nrow = 1,ncol = 8,
          labels=c("I","J",'K',"L","M",'N',"O","P"))
dev.off()

ggarrange(p_TU_fig1g,p_TUY_pearson,ncol=2,nrow=1,labels=c("D","E"),widths=c(0.6,0.8), heights=c(0.5,0.5),align ="hv")

#Fig S1 H) GFP小提琴图------------
#TGY-----
TGY_pearson_nowt <- read.table("TGY_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(TGY_pearson_nowt) <- c("x","y","pearson")

TGY_pearson_nowt$x_sample <- TGY_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TGY_pearson_nowt$x_time <- TGY_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
TGY_pearson_nowt$y_sample <- TGY_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TGY_pearson_nowt$y_time <- TGY_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

TGY_pearson_nowt_group <- lapply(1:nrow(TGY_pearson_nowt), function(x){
  data <- TGY_pearson_nowt[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}

}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'TGY')


#TGS-----
TGS_pearson_nowt <- read.table("TGS_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(TGS_pearson_nowt) <- c("x","y","pearson")

TGS_pearson_nowt$x_sample <- TGS_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TGS_pearson_nowt$x_time <- TGS_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
TGS_pearson_nowt$y_sample <- TGS_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TGS_pearson_nowt$y_time <- TGS_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

TGS_pearson_nowt_group <- lapply(1:nrow(TGS_pearson_nowt), function(x){
  data <- TGS_pearson_nowt[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'TGS')

#AGY-----
AGY_pearson_nowt <- read.table("AGY_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(AGY_pearson_nowt) <- c("x","y","pearson")

AGY_pearson_nowt$x_sample <- AGY_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AGY_pearson_nowt$x_time <- AGY_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
AGY_pearson_nowt$y_sample <- AGY_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AGY_pearson_nowt$y_time <- AGY_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

AGY_pearson_nowt_group <- lapply(1:nrow(AGY_pearson_nowt), function(x){
  data <- AGY_pearson_nowt[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'AGY')

#AGS-----
AGS_pearson_nowt <- read.table("AGS_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(AGS_pearson_nowt) <- c("x","y","pearson")

AGS_pearson_nowt$x_sample <- AGS_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AGS_pearson_nowt$x_time <- AGS_pearson_nowt$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
AGS_pearson_nowt$y_sample <- AGS_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AGS_pearson_nowt$y_time <- AGS_pearson_nowt$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

AGS_pearson_nowt_group <- lapply(1:nrow(AGS_pearson_nowt), function(x){
  data <- AGS_pearson_nowt[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'AGS')


# 四个样本合并画图-----

TA_G_YS_pearson_nowt <- rbind(TGY_pearson_nowt_group,TGS_pearson_nowt_group,AGY_pearson_nowt_group,AGS_pearson_nowt_group)


library(RColorBrewer)
my_orange = brewer.pal(n = 8, "Set2")[3:6]

TA_G_YS_pearson_nowt %>%
  mutate(x = factor(group,levels = c('biol replicates','diff time of the same sample','diff time of the diff sample'),order=T))%>%
  mutate(set = factor(sam,levels = c('TGY','TGS','AGY','AGS'),order=T))%>%
  ggplot(aes(x = x, y = pearson, fill = set)) +
  scale_fill_manual(values = my_orange)+
  geom_violin(alpha=0.8,width=0.6,scale = 'width') + ylab('pearson') + xlab('') +
  #facet_grid(.~set, scales = "free", switch = "x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  style.print()

#facet_grid(.~set, scales = "free", switch = "x", space = "free_x") +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# 11-20 Fig S1H 用TG数据-----------

TA_G_YS_pearson <- rbind(TGY_pearson_nowt_group,AGY_pearson_nowt_group)

TA_G_YS_pearson<- TA_G_YS_pearson %>%
  mutate(new_group = ifelse(group == 'biol replicates','biol repeat','not biol repeat'))

TA_G_YS_pearson %>%
  ggplot(aes(x=new_group, y=pearson, fill = new_group)) +
  geom_boxplot(alpha=0.8,width=0.6)+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  style.print()



my_color2 = brewer.pal(n = 12, "Paired")[2:3]
ggplot(TA_G_YS_pearson, aes(x=pearson, fill=new_group)) +
  geom_histogram(data = filter(TA_G_YS_pearson,new_group == 'biol repeat'), aes(y=..count../sum(..count..)),bins = 12, position = 'identity',alpha=0.5) +
  geom_histogram(data = filter(TA_G_YS_pearson,new_group == 'not biol repeat'), aes(y=..count../sum(..count..)),bins = 12, position = 'identity',alpha=0.5) +
  #scale_fill_manual(values= my_color2)+
  #stat_bin(bins = 15,position = 'dodge')+
  #stat_bin(bins = 10,position="dodge") +
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


TA_G_YS_pearson_bio <- filter(TA_G_YS_pearson,new_group == 'biol repeat')
TA_G_YS_pearson_nbio <- filter(TA_G_YS_pearson,new_group == 'not biol repeat')

filter(TA_G_YS_pearson_nbio,pearson>0.78&pearson<=0.8) %>% nrow()

test3 <- lapply(1:12, function(x){
  num<-0.76+(x-1)*0.02
  num1<- 0.76+x*0.02
  tet <- filter(TA_G_YS_pearson_nbio,pearson>num&pearson<=num1) %>% nrow()
  data <- data.frame(median=(num+num1)/2,freq=tet/216,group='Not biological repeat')
}) %>%rbind.fill()

test4 <- lapply(1:12, function(x){
  num<-0.76+(x-1)*0.02
  num1<- 0.76+x*0.02
  tet <- filter(TA_G_YS_pearson_bio,pearson>num&pearson<=num1) %>% nrow()
  data <- data.frame(median=(num+num1)/2,freq=tet/48,group='Biological repeat')
}) %>%rbind.fill()

test5 <- rbind(test3,test4)

p_figs1h <- s

test5 %>%
  mutate(x=factor(median,levels = c('0.77','0.79','0.81','0.83','0.85','0.87','0.89','0.91','0.93','0.95','0.97','0.99'),ordered = T)) %>%
  ggplot(aes(x=x,y=freq,group=group,fill=group)) +
  geom_bar(stat="identity",position=position_dodge(width=0.5,preserve="single")) +
  theme_bw() + xlab('Pearson correlation') +ylab('Frequency') +
  style.print()+ scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  scale_x_discrete(labels=c('(0.76,0.78]','(0.78,0.8]','(0.8,0.82]','(0.82,0.84]','(0.84,0.86]','(0.86,0.88]','(0.88,0.9]','(0.9,0.92]',
                            '(0.92,0.94]','(0.94,0.96]','(0.96,0.98]','(0.98,1]'))+
  theme(axis.text=element_text(size=18),
        axis.title.x =element_text(size=22),
        axis.title.y=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22),axis.text.x = element_text(angle = 45, hjust = 1))

#--------------------
#测试分6个bar   more exquisite!!! just in my view

test6 <- lapply(1:4, function(x){
    num<-0.88+(x-1)*0.03
    num1<- 0.88+x*0.03
    tet <- filter(TA_G_YS_pearson_nbio,pearson>num&pearson<=num1) %>% nrow()
    data <- data.frame(median=(num+num1)/2,freq=tet/216,group='Not biological repeat')
  }) %>%rbind.fill()

test7 <- lapply(1:4, function(x){
  num<-0.88+(x-1)*0.03
  num1<- 0.88+x*0.03
  tet <- filter(TA_G_YS_pearson_bio,pearson>num&pearson<=num1) %>% nrow()
  data <- data.frame(median=(num+num1)/2,freq=tet/48,group='Biological repeat')
}) %>%rbind.fill()


test9 <- lapply(1:5, function(x){
  num<-0.88+(x-1)*0.03
  num1<- 0.88+x*0.03
  data <- data.frame(median=num,freq=0,group='Biological repeat')
}) %>%rbind.fill()

test10 <- lapply(1:5, function(x){
  num<-0.88+(x-1)*0.03
  num1<- 0.88+x*0.03
  data <- data.frame(median=num,freq=0,group='Not biological repeat')
}) %>%rbind.fill()


test8 <- rbind(test6,test7,test9,test10)


p_figs1h <- test8 %>%
  mutate(x=factor(median,levels = c('0.88','0.895','0.91','0.925','0.94',
                                    '0.955','0.97','0.985','1'),ordered = T)) %>%
  ggplot(aes(x=x,y=freq,group=group,fill=group)) +
  geom_bar(stat="identity",position=position_dodge(width=0.65,preserve="single"),
           width = 1.3) +
  theme_bw() + xlab('Pearson correlation') +ylab('Frequency') +
  style.print()+ scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  scale_x_discrete(breaks= c('0.88','0.91','0.94','0.97','1'),
                   labels=c('0.88','0.91','0.94','0.97','1'))+
  scale_fill_manual(labels=c('Different time','Same time'),values = c('#00BFC4','#F8766D'))+
  labs(fill='Samples between')+
  theme(axis.text=element_text(size=18),
        axis.title.x =element_text(size=22),
        axis.title.y=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22))

pdf("~/figS1h.pdf",width = 14,height = 10)
p_figs1h
dev.off()



#--------------------------------

test5 %>%
  #mutate(x=factor(median,levels = c('0.79','0.81','0.83','0.85','0.87','0.89','0.91','0.93','0.95','0.97','0.99'),ordered = T)) %>%
  ggplot(aes(x=factor(median),y=freq,group=group,fill=group)) +
  geom_bar(stat="identity",position=position_dodge(width=0.5)) +
  theme_bw() + xlab('Pearson correlation') +ylab('Frequency') +
  style.print()+ scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  scale_x_discrete(breaks=seq(0,11,0.5),
                   labels=c('0.78','0.8','0.82','0.84','0.86','0.88','0.9','0.92',
                            '0.94','0.78','0.8','0.82','0.84','0.86','0.88','0.9','0.92',
                            '0.94','c','c','d','f','s'))+
  theme(axis.text=element_text(size=20),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        legend.title = element_text(size=24))

pdf("~/figS1h.pdf",width = 14,height = 10)
p_figs1h
dev.off()


ggplot(TA_G_YS_pearson_bio,aes(x=pearson)) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 12, position = 'identity',alpha=0.5) +
  stat_bin(bins = 12, geom="text", aes(label=..count../sum(..count..)), vjust=-1.5) +
  #scale_fill_manual(values= my_color2)+
  #stat_bin(bins = 15,position = 'dodge')+
  #stat_bin(bins = 10,position="dodge") +
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


my_color2 = brewer.pal(n = 12, "Paired")[2:3]
ggplot(TA_G_YS_pearson) +
  geom_histogram(aes(x=pearson, fill=new_group,y=..count../sum(..count..)),bins = 12, position = 'dodge',alpha=0.5) +
  #scale_fill_manual(values= my_color2)+
  #stat_bin(bins = 15,position = 'dodge')+
  #stat_bin(bins = 10,position="dodge") +
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5)) + ylim(0,0.6)















#Fig S1 H) URA3小提琴图------------
#TUY-----
TUY_pearson <- read.csv("TUY_pearson.csv",header = F,stringsAsFactors = F)
names(TUY_pearson) <- c("x","y","pearson")

TUY_pearson$x_sample <- TUY_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TUY_pearson$x_time <- TUY_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
TUY_pearson$y_sample <- TUY_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TUY_pearson$y_time <- TUY_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

TUY_pearson_group <- lapply(1:nrow(TUY_pearson), function(x){
  data <- TUY_pearson[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'TUY')


#TUS-----
TUS_pearson <- read.csv("TUS_pearson.csv",header = F, stringsAsFactors = F)
names(TUS_pearson) <- c("x","y","pearson")

TUS_pearson$x_sample <- TUS_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TUS_pearson$x_time <- TUS_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
TUS_pearson$y_sample <- TUS_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
TUS_pearson$y_time <- TUS_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

TUS_pearson_group <- lapply(1:nrow(TUS_pearson), function(x){
  data <- TUS_pearson[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'TUS')

#AUY-----
AUY_pearson <- read.csv("AUY_pearson.csv",header = F, stringsAsFactors = F)
names(AUY_pearson) <- c("x","y","pearson")

AUY_pearson$x_sample <- AUY_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AUY_pearson$x_time <- AUY_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
AUY_pearson$y_sample <- AUY_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AUY_pearson$y_time <- AUY_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

AUY_pearson_group <- lapply(1:nrow(AUY_pearson), function(x){
  data <- AUY_pearson[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'AUY')

#AUS-----
AUS_pearson <- read.csv("AUS_pearson.csv",header = F, stringsAsFactors = F)
names(AUS_pearson) <- c("x","y","pearson")

AUS_pearson$x_sample <- AUS_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AUS_pearson$x_time <- AUS_pearson$x %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()
AUS_pearson$y_sample <- AUS_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[1]}) %>%
  unlist()
AUS_pearson$y_time <- AUS_pearson$y %>%
  strsplit("-") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

AUS_pearson_group <- lapply(1:nrow(AUS_pearson), function(x){
  data <- AUS_pearson[x,]
  if (data$x_time == data$y_time & data$x_sample != data$y_sample) {
    dat <- data %>% mutate(group='biol replicates')
  }else if(data$x_time != data$y_time & data$x_sample == data$y_sample){
    dat <- data %>% mutate(group='diff time of the same sample')
  }else if(data$x_time != data$y_time & data$x_sample != data$y_sample){
    dat <- data %>% mutate(group='diff time of the diff sample')
  }else{dat <- data %>% mutate(group='0')}
  
}) %>% rbind.fill() %>% filter(group!="0") %>% mutate(sam = 'AUS')

AUS_pearson_group <- dplyr::filter(AUS_pearson_group,AUS_pearson_group$x!='s3-d7' & AUS_pearson_group$y!='s3-d7')


# 四个样本合并画图-----

TA_U_YS_pearson <- rbind(TUY_pearson_group,TUS_pearson_group,AUY_pearson_group,AUS_pearson_group)


library(RColorBrewer)
my_orange = brewer.pal(n = 8, "Set2")[3:6]

TA_U_YS_pearson %>%
  mutate(x = factor(group,levels = c('biol replicates','diff time of the same sample','diff time of the diff sample'),order=T))%>%
  mutate(set = factor(sam,levels = c('TUY','TUS','AUY','AUS'),order=T))%>%
  ggplot(aes(x = x, y = pearson,fill=set)) +
  scale_fill_manual(values = my_orange)+
  geom_point(aes(group = set)) + ylab('pearson') + xlab('') +
  geom_boxplot(alpha=0.8,width=0.6,scale = 'width')+
  #facet_grid(.~set, scales = "free", switch = "x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  style.print()

TA_U_YS_pearson<- TA_U_YS_pearson %>%
  mutate(new_group = ifelse(group == 'biol replicates','biol repeat','not biol repeat'))

TA_U_YS_pearson %>%
  ggplot(aes(x=new_group, y=pearson, fill = new_group)) +
  geom_boxplot(alpha=0.8,width=0.6)+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  style.print()



my_color2 = brewer.pal(n = 12, "Paired")[2:3]
ggplot(TA_U_YS_pearson, aes(x=pearson, fill=new_group)) +
geom_histogram(data = filter(TA_U_YS_pearson,new_group == 'biol repeat'), aes(y=..count../sum(..count..)),bins = 12, position = 'identity',alpha=0.5) +
  geom_histogram(data = filter(TA_U_YS_pearson,new_group == 'not biol repeat'), aes(y=..count../sum(..count..)),bins = 12, position = 'identity',alpha=0.5) +
  #scale_fill_manual(values= my_color2)+
  #stat_bin(bins = 15,position = 'dodge')+
  #stat_bin(bins = 10,position="dodge") +
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


#1113 new S1H 叠0.5的bar图 -------------

TA_U_YS_pearson_bio <- filter(TA_U_YS_pearson,new_group == 'biol repeat')
TA_U_YS_pearson_nbio <- filter(TA_U_YS_pearson,new_group == 'not biol repeat')

filter(TA_U_YS_pearson_nbio,pearson>0.78&pearson<=0.8) %>% nrow()

test <- lapply(1:11, function(x){
  num<-0.78+(x-1)*0.02
  num1<- 0.78+x*0.02
  tet <- filter(TA_U_YS_pearson_nbio,pearson>num&pearson<=num1) %>% nrow()
  data <- data.frame(median=(num+num1)/2,freq=tet/414,group='Not biological repeat')
}) %>%rbind.fill()

test1 <- lapply(1:11, function(x){
  num<-0.78+(x-1)*0.02
  num1<- 0.78+x*0.02
  tet <- filter(TA_U_YS_pearson_bio,pearson>num&pearson<=num1) %>% nrow()
  data <- data.frame(median=(num+num1)/2,freq=tet/92,group='Biological repeat')
}) %>%rbind.fill()

test2 <- rbind(test,test1)

p_figs1h <- test2 %>%
  mutate(x=factor(median,levels = c('0.79','0.81','0.83','0.85','0.87','0.89','0.91','0.93','0.95','0.97','0.99'),ordered = T)) %>%
  ggplot(aes(x=x,y=freq,group=group,fill=group)) +
  geom_bar(stat="identity",position=position_dodge(width=0.5,preserve="single")) +
  theme_bw() + xlab('Pearson correlation') +ylab('Frequency') +
  style.print()+ scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  scale_x_discrete(labels=c('(0.78,0.8]','(0.8,0.82]','(0.82,0.84]','(0.84,0.86]','(0.86,0.88]','(0.88,0.9]','(0.9,0.92]',
                                                      '(0.92,0.94]','(0.94,0.96]','(0.96,0.98]','(0.98,1]'))+
  theme(axis.text=element_text(size=18),
        axis.title.x =element_text(size=22),
        axis.title.y=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22),axis.text.x = element_text(angle = 45, hjust = 1))


test2 %>%
  #mutate(x=factor(median,levels = c('0.79','0.81','0.83','0.85','0.87','0.89','0.91','0.93','0.95','0.97','0.99'),ordered = T)) %>%
  ggplot(aes(x=factor(median),y=freq,group=group,fill=group)) +
  geom_bar(stat="identity",position=position_dodge(width=0.5)) +
  theme_bw() + xlab('Pearson correlation') +ylab('Frequency') +
  style.print()+ scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  scale_x_discrete(breaks=seq(0,11,0.5),
                   labels=c('0.78','0.8','0.82','0.84','0.86','0.88','0.9','0.92',
                            '0.94','0.78','0.8','0.82','0.84','0.86','0.88','0.9','0.92',
                            '0.94','c','c','d','f','s'))+
  theme(axis.text=element_text(size=20),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        legend.title = element_text(size=24))

pdf("~/figS1h.pdf",width = 14,height = 10)
p_figs1h
dev.off()


ggplot(TA_U_YS_pearson_bio,aes(x=pearson)) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 12, position = 'identity',alpha=0.5) +
  stat_bin(bins = 12, geom="text", aes(label=..count../sum(..count..)), vjust=-1.5) +
  #scale_fill_manual(values= my_color2)+
  #stat_bin(bins = 15,position = 'dodge')+
  #stat_bin(bins = 10,position="dodge") +
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


my_color2 = brewer.pal(n = 12, "Paired")[2:3]
ggplot(TA_U_YS_pearson) +
  geom_histogram(aes(x=pearson, fill=new_group,y=..count../sum(..count..)),bins = 12, position = 'dodge',alpha=0.5) +
  #scale_fill_manual(values= my_color2)+
  #stat_bin(bins = 15,position = 'dodge')+
  #stat_bin(bins = 10,position="dodge") +
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5)) + ylim(0,0.6)




# S1H 累计分布曲线

ggplot(TA_U_YS_pearson, aes(x=pearson)) +
 stat_ecdf(aes(group=new_group,colour = new_group)) +  #这里也可以画geom_line stat='ecdf'即可
  theme_bw() + xlab('Pearson') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))




#FigS2 CD)------
#TGY-----
TGY_day7_plotfil <- TGY_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried')) %>%
  mutate(eff_size = fitness_TGY-1)%>%
  mutate(abs_eff_size = abs(eff_size))

my_color = brewer.pal(n = 8, "Set2")[1:2]


#累积频率分布直方图（同义非同义）-------
#p_syno<- 
  
figs2c_TGY <-ggplot(TGY_day7_plotfil, aes(x=eff_size, colour=group)) +
  scale_fill_manual(values= my_color)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('s') +ylab('Cumulative frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5),
        legend.position = c(0.9,0.1),
        legend.justification  = c(0.9,0.1))+
  guides(colour=guide_legend(title=NULL))+
  #labs(title = 'TGY synon-nonsynon')+
  annotate('text',parse = TRUE, label = " ~ italic(P)< ~10^-18", x = 0.01, y = 0.75,size=10)

wilcox.test(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous')],
            TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],
            alternative = 'greater')$p.value

# TGY qq-----
s_nonsy_TGY <- sort(quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_TGY <- sort(quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_TGY_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_TGY, y=s_sy_TGY),size=5) +
  style.print() +
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  " ~ italic(P)< ~10^-15", x = -0.005, y = 0.017,size=50)


wilcox.test(TGY_sy_box$eff_size[TGY_sy_box$group=='synonymous'],TGY_sy_box$eff_size[TGY_sy_box$group=='nonsynonymous'],alternative = 'greater')$p.value
1.381605e-16


#__________________________________暂时不用（分了百分位后的boxplot）
TGY_sy <- data.frame(s_sy)%>%
  mutate(type = 'sy')
names(TGY_sy) = c('eff_s','type')
TGY_nonsy <- data.frame(s_nonsy)%>%
  mutate(tyoe = 'nonsy')
names(TGY_nonsy) = c('eff_s','type')
TGY_sy_nonsy <- rbind(TGY_sy,TGY_nonsy) %>%
  mutate(group = 'TGY')

ggplot(TGY_sy_nonsy) + 
  geom_boxplot(aes(x=type, y=eff_s)) +
  #geom_abline(intercept = 0, slope = 1 ) +
  #geom_hline(aes(yintercept=0), linetype = "dashed") +
  #geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  labs(title = 'TGY synon-nonsynon')
#__________________________________






s_nonsy_abs <- sort(quantile(TGY_day7_plotfil$abs_eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_abs <- sort(quantile(TGY_day7_plotfil$abs_eff_size[which(TGY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))
lenx <- length(s_nonsy_abs)
leny <- length(s_sy_abs)
if (leny < lenx)s_nonsy <- approx(1L:lenx, s_nonsy_abs, n = leny)$y
if (leny > lenx)s_sy <- approx(1L:leny, s_sy_abs, n = lenx)$y
p_TGY_sy_qq_abs<- ggplot() + 
  geom_point(aes(x=s_nonsy_abs, y=s_sy_abs)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  labs(title = 'TGY synon-nonsynon abs')

ggarrange(p_TGY_sy_qq,p_TGY_sy_qq_abs,ncol = 2,nrow = 1)


wilcox.test(TGY_day7_plotfil$abs_eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],TGY_day7_plotfil$abs_eff_size[which(TGY_day7_plotfil$group=='synonymous')],alternative = 'less')
wilcox.test(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group == 'synonymous')] )


quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous')])
quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')])
wilcox.test(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous' & TGY_day7_plotfil$eff_size >quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],0.95))],
            TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous'& TGY_day7_plotfil$eff_size >quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous')],0.95))],alternative ="greater")

wilcox.test(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous' & TGY_day7_plotfil$eff_size <quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='nonsynonymous')],0.05))],
            TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous'& TGY_day7_plotfil$eff_size <quantile(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$group=='synonymous')],0.05))],alternative ="less")


my_color2 = brewer.pal(n = 12, "Paired")[2:3]

# 累积频率分布直方图（内部外部）
ggplot(TGY_day7_plotfil, aes(x=eff_size, colour=type)) +
  scale_fill_manual(values= my_color2)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


# histogram（内部外部）
ggplot(TGY_day7_plotfil, aes(x=eff_size, fill=type)) +
  scale_fill_manual(values= my_color2)+
  stat_bin(alpha=0.8, position="identity",bins = 40) +
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

# histogram只要轮廓（内部外部）
ggplot(TGY_day7_plotfil, aes(x=eff_size, colour=type)) +
  stat_bin(aes(y=..count../sum(..count..)),bins = 40,geom="step",position = 'identity') +
  theme_bw() + xlab('eff_size') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

ggarrange(p_syno,p_surcore,ncol = 2,nrow = 1)

wilcox.test(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$type=='surfaced')],TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$type=='buried')])
median(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$type=='surfaced')])
median(TGY_day7_plotfil$eff_size[which(TGY_day7_plotfil$type=='buried')])


#AGY-----
AGY_day7_plotfil <- AGY_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_AGY-1)%>%
  mutate(abs_eff_size = abs(eff_size))

my_color = brewer.pal(n = 8, "Set2")[1:2]

#累积频率分布直方图（同义非同义）------
#p_syno<- 

figs2c_AGY <-ggplot(AGY_day7_plotfil, aes(x=eff_size, colour=group)) +
  scale_fill_manual(values= my_color)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Cumulative frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5),
        legend.position = c(0.9,0.1),
        legend.justification  = c(0.9,0.1))+
  guides(colour=guide_legend(title=NULL))+
  labs(title = 'AGY synon-nonsynon')+
  annotate('text', label = 'Wilcox test(abs) P = 2.63e-06(nonsynon)', x = 0.01, y = 0.75,size=6)

#AGY qq ------
s_nonsy_AGY <- sort(quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_AGY <- sort(quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_AGY_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_AGY, y=s_sy_AGY),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  " ~ italic(P) == 0.215", x = -0.005, y = 0.017,size=50)

wilcox.test(AGY_sy_box$eff_size[AGY_sy_box$group=='synonymous'],AGY_sy_box$eff_size[AGY_sy_box$group=='nonsynonymous'],alternative = 'greater')$p.value
0.2146542



s_nonsy_abs <- sort(quantile(AGY_day7_plotfil$abs_eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_abs <- sort(quantile(AGY_day7_plotfil$abs_eff_size[which(AGY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))
lenx <- length(s_nonsy_abs)
leny <- length(s_sy_abs)
if (leny < lenx)s_nonsy <- approx(1L:lenx, s_nonsy_abs, n = leny)$y
if (leny > lenx)s_sy <- approx(1L:leny, s_sy_abs, n = lenx)$y
p_AGY_sy_qq_abs<- ggplot() + 
  geom_point(aes(x=s_nonsy_abs, y=s_sy_abs)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  labs(title = 'AGY synon-nonsynon abs')

ggarrange(p_AGY_sy_qq,p_AGY_sy_qq_abs,ncol = 2,nrow = 1)




wilcox.test(AGY_day7_plotfil$abs_eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')],AGY_day7_plotfil$abs_eff_size[which(AGY_day7_plotfil$group=='synonymous')],alternative ="greater")
wilcox.test(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')],AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group == 'synonymous')],alternative ="less")


quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='synonymous')])
quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')])
wilcox.test(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous' & AGY_day7_plotfil$eff_size >quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')],0.9))],
            AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='synonymous'& AGY_day7_plotfil$eff_size >quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='synonymous')],0.9))],alternative ="greater")

wilcox.test(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous' & AGY_day7_plotfil$eff_size <quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='nonsynonymous')],0.05))],
            AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='synonymous'& AGY_day7_plotfil$eff_size <quantile(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$group=='synonymous')],0.05))],alternative ="less")


my_color2 = brewer.pal(n = 12, "Paired")[2:3]

# 累积频率分布直方图（内部外部）
ggplot(AGY_day7_plotfil, aes(x=eff_size, colour=type)) +
  scale_fill_manual(values= my_color2)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


# histogram（内部外部）
ggplot(AGY_day7_plotfil, aes(x=eff_size, fill=type)) +
  scale_fill_manual(values= my_color2)+
  stat_bin(alpha=0.8, position="identity",bins = 40) +
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

# histogram只要轮廓（内部外部）
ggplot(AGY_day7_plotfil, aes(x=eff_size, colour=type)) +
  stat_bin(aes(y=..count../sum(..count..)),bins = 40,geom="step",position = 'identity') +
  theme_bw() + xlab('eff_size') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

ggarrange(p_syno,p_surcore,ncol = 2,nrow = 1)

wilcox.test(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$type=='surfaced')],AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$type=='buried')])
median(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$type=='surfaced')])
median(AGY_day7_plotfil$eff_size[which(AGY_day7_plotfil$type=='buried')])

#TGS-----
TGS_day7_plotfil <- TGS_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried')) %>%
  mutate(eff_size = fitness_TGS-1)%>%
  mutate(abs_eff_size = abs(eff_size))

#TGS qq ------
s_nonsy_TGS <- sort(quantile(TGS_day7_plotfil$eff_size[which(TGS_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_TGS <- sort(quantile(TGS_day7_plotfil$eff_size[which(TGS_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_TGS_sy_qq<- ggplot() + 
  geom_point(aes(x=s_nonsy_TGS, y=s_sy_TGS)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed") +
  style.print()+
  scale_x_continuous(limits = c(-0.03,0.03),breaks = c(-0.03,0,0.03))+
  scale_y_continuous(limits = c(-0.03,0.03),breaks = c(-0.03,0,0.03))+
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=70,hjust = 0.8),axis.text.y = element_text(size=70),
        axis.title = element_text(size=70))+
  annotate('text',parse = TRUE,label =  " ~ italic(P)< ~10^-14", x = -0.015, y = 0.03,size=25)


wilcox.test(TGS_sy_box$eff_size[TGS_sy_box$group=='synonymous'],TGS_sy_box$eff_size[TGS_sy_box$group=='nonsynonymous'],alternative = 'greater')$p.value
4.885e-15

#AGS-----
AGS_day7_plotfil <- AGS_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried')) %>%
  mutate(eff_size = fitness_AGS-1)%>%
  mutate(abs_eff_size = abs(eff_size))

#AGS qq ------
s_nonsy_AGS <- sort(quantile(AGS_day7_plotfil$eff_size[which(AGS_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_AGS <- sort(quantile(AGS_day7_plotfil$eff_size[which(AGS_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_AGS_sy_qq<- ggplot() + 
  geom_point(aes(x=s_nonsy_AGS, y=s_sy_AGS)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  scale_x_continuous(limits = c(-0.05,0.05),breaks = c(-0.05,0,0.05))+
  scale_y_continuous(limits = c(-0.05,0.05),breaks = c(-0.05,0,0.05))+
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=70,hjust = 0.8),axis.text.y = element_text(size=70),
        axis.title = element_text(size=70))+
  annotate('text',parse = TRUE,label =  " ~ italic(P) == 0.064", x = -0.025, y = 0.05,size=25)


wilcox.test(AGS_sy_box$eff_size[AGS_sy_box$group=='synonymous'],AGS_sy_box$eff_size[AGS_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.06373



#TUY-----
TUY_day7_plotfil <- TUY_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_TUY-1) %>%
  mutate(abs_eff_size = abs(eff_size))


#累积频率分布直方图（同义非同义）------
my_color = brewer.pal(n = 8, "Set2")[1:2]

figs2c_TUY <-ggplot(TUY_day7_plotfil, aes(x=eff_size, colour=group)) +
  scale_fill_manual(values= my_color)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Cumulative frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5),
        legend.position = c(0.9,0.1),
        legend.justification  = c(0.9,0.1))+
  guides(colour=guide_legend(title=NULL))+
  labs(title = 'TUY synon-nonsynon')+
  annotate('text', label = 'Wilcox test(abs) P = 0.0045(nonsynon)', x = 0.02, y = 0.75,size=6)


# TUY qq-----
s_nonsy_TUY <- sort(quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_TUY <- sort(quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_TUY_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_TUY, y=s_sy_TUY),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  "~ italic(P)< ~0.02", x = -0.005, y = 0.017,size=50)


wilcox.test(TUY_sy_box$eff_size[TUY_sy_box$group=='synonymous'],TUY_sy_box$eff_size[TUY_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.01001


s_nonsy_abs <- sort(quantile(TUY_day7_plotfil$abs_eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_abs <- sort(quantile(TUY_day7_plotfil$abs_eff_size[which(TUY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))
lenx <- length(s_nonsy_abs)
leny <- length(s_sy_abs)
if (leny < lenx)s_nonsy <- approx(1L:lenx, s_nonsy_abs, n = leny)$y
if (leny > lenx)s_sy <- approx(1L:leny, s_sy_abs, n = lenx)$y
p_TUY_sy_qq_abs<- ggplot() + 
  geom_point(aes(x=s_nonsy_abs, y=s_sy_abs)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  labs(title = 'TUY synon-nonsynon abs')

ggarrange(p_TUY_sy_qq,p_TUY_sy_qq_abs,ncol = 2,nrow = 1)

wilcox.test(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')],TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous')],alternative ="less")
wilcox.test(TUY_day7_plotfil$abs_eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')],TUY_day7_plotfil$abs_eff_size[which(TUY_day7_plotfil$group=='synonymous')],alternative ="greater")

#quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous')])
#quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')])
wilcox.test(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous' & TUY_day7_plotfil$eff_size >quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')],0.9))],
            TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous'& TUY_day7_plotfil$eff_size >quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous')],0.9))],alternative ="greater")

wilcox.test(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous' & TUY_day7_plotfil$eff_size <quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='nonsynonymous')],0.1))],
            TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous'& TUY_day7_plotfil$eff_size <quantile(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$group=='synonymous')],0.1))],alternative ="less")


my_color2 = brewer.pal(n = 12, "Paired")[2:3]
# 累积频率分布直方图（内部外部）
ggplot(TUY_day7_plotfil, aes(x=eff_size, colour=type)) +
  scale_fill_manual(values= my_color2)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

# histogram（内部外部）
ggplot(TUY_day7_plotfil, aes(x=eff_size, fill=type)) +
  scale_fill_manual(values= my_color2)+
  stat_bin(alpha=0.8, position="identity",bins = 40) +
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

# histogram只要轮廓（内部外部）
ggplot(TUY_day7_plotfil, aes(x=eff_size, colour=type)) +
  stat_bin(aes(y=..count../sum(..count..)),bins = 40,geom="step",position = 'identity') +
  theme_bw() + xlab('eff_size') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

wilcox.test(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$type=='surfaced')],TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$type=='buried')],alternative ="greater" )
median(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$type=='surfaced')])
median(TUY_day7_plotfil$eff_size[which(TUY_day7_plotfil$type=='buried')])

ggarrange(p_syno,p_surcore,ncol = 2,nrow = 1)

#AUY-----
AUY_day7_plotfil <- AUY_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_AUY-1)%>%
  mutate(abs_eff_size = abs(eff_size))

#累积频率分布直方图（同义非同义）------
my_color = brewer.pal(n = 8, "Set2")[1:2]

figs2c_AUY <-ggplot(AUY_day7_plotfil, aes(x=eff_size, colour=group)) +
  scale_fill_manual(values= my_color)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Cumulative frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5),
        legend.position = c(0.9,0.1),
        legend.justification  = c(0.9,0.1))+
  guides(colour=guide_legend(title=NULL))+
  labs(title = 'AUY synon-nonsynon')+
  annotate('text', label = 'Wilcox test(abs) P = 0.0027(synon)', x = 0.03, y = 0.75,size=6)

#AUY qq -----
s_nonsy_AUY <- sort(quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_AUY <- sort(quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_AUY_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_AUY, y=s_sy_AUY),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  "~ italic(P)< ~0.02", x = -0.005, y = 0.017,size=50)


wilcox.test(AUY_sy_box$eff_size[AUY_sy_box$group=='synonymous'],AUY_sy_box$eff_size[AUY_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.01837

wilcox.test(AUY_day7_plotfil$abs_eff_size[which(AUY_day7_plotfil$group=='nonsynonymous')],AUY_day7_plotfil$abs_eff_size[which(AUY_day7_plotfil$group=='synonymous')],alternative ="less")
wilcox.test(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous')],AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group == 'synonymous')],alternative ="less" )


quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='synonymous')])
quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous')])
wilcox.test(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous' & AUY_day7_plotfil$eff_size >quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous')],0.9))],
            AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='synonymous'& AUY_day7_plotfil$eff_size >quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='synonymous')],0.9))],alternative ="greater")

wilcox.test(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous' & AUY_day7_plotfil$eff_size <quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='nonsynonymous')],0.05))],
            AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='synonymous'& AUY_day7_plotfil$eff_size <quantile(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$group=='synonymous')],0.05))],alternative ="less")


ggarrange(figs2c_TGY,figs2c_AGY,figs2c_TUY,figs2c_AUY,ncol = 2, nrow = 2)



my_color2 = brewer.pal(n = 12, "Paired")[2:3]

# 累积频率分布直方图（内部外部）
ggplot(AUY_day7_plotfil, aes(x=eff_size, colour=type)) +
  scale_fill_manual(values= my_color2)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))


# histogram（内部外部）
ggplot(AUY_day7_plotfil, aes(x=eff_size, fill=type)) +
  scale_fill_manual(values= my_color2)+
  stat_bin(alpha=0.8, position="identity",bins = 40) +
  theme_bw() + xlab('eff_size') +ylab('Count') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

# histogram只要轮廓（内部外部）
ggplot(AUY_day7_plotfil, aes(x=eff_size, colour=type)) +
  stat_bin(aes(y=..count../sum(..count..)),bins = 40,geom="step",position = 'identity') +
  theme_bw() + xlab('eff_size') +ylab('Frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5))

ggarrange(p_surcore,p_surcore2,ncol = 2,nrow = 1)

wilcox.test(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$type=='surfaced')],AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$type=='buried')])
median(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$type=='surfaced')])
median(AUY_day7_plotfil$eff_size[which(AUY_day7_plotfil$type=='buried')])


#TUS-----
TUS_day7_plotfil <- TUS_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_TUS-1) %>%
  mutate(abs_eff_size = abs(eff_size))

#TUS qq -----
s_nonsy_TUS <- sort(quantile(TUS_day7_plotfil$eff_size[which(TUS_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_TUS <- sort(quantile(TUS_day7_plotfil$eff_size[which(TUS_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_TUS_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_TUS, y=s_sy_TUS),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  "~ italic(P)< ~10^-18", x = -0.005, y = 0.017,size=50)



wilcox.test(TUS_sy_box$eff_size[TUS_sy_box$group=='synonymous'],TUS_sy_box$eff_size[TUS_sy_box$group=='nonsynonymous'],alternative = 'greater')$p.value
2.02751e-19


#AUS-----
AUS_day7_plotfil <- AUS_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_AUS-1)%>%
  mutate(abs_eff_size = abs(eff_size))

#AUS qq -----
s_nonsy_AUS <- sort(quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_AUS <- sort(quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_AUS_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_AUS, y=s_sy_AUS),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  "~ italic(P)< ~10^-11", x = -0.005, y = 0.017,size=50)



wilcox.test(AUS_sy_box$eff_size[AUS_sy_box$group=='synonymous'],AUS_sy_box$eff_size[AUS_sy_box$group=='nonsynonymous'],alternative = 'greater')
2.383477e-12


#累积频率分布直方图（同义非同义）------
my_color = brewer.pal(n = 8, "Set2")[1:2]

figs2c_AUS <-ggplot(AUS_day7_plotfil, aes(x=eff_size, colour=group)) +
  scale_fill_manual(values= my_color)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Cumulative frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5),
        legend.position = c(0.9,0.1),
        legend.justification  = c(0.9,0.1))+
  guides(colour=guide_legend(title=NULL))+
  labs(title = 'AUS synon-nonsynon')
#annotate('text', label = 'Wilcox test(abs) P = 0.0027(synon)', x = 0.03, y = 0.75,size=6)

s_nonsy <- sort(quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy <- sort(quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))
lenx <- length(s_nonsy)
leny <- length(s_sy)
if (leny < lenx)s_nonsy <- approx(1L:lenx, s_nonsy, n = leny)$y
if (leny > lenx)s_sy <- approx(1L:leny, s_sy, n = lenx)$y

p_AUS_sy_qq<- ggplot() + 
  geom_point(aes(x=s_nonsy, y=s_sy)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  labs(title = 'AUS synon-nonsynon')

s_nonsy_abs <- sort(quantile(AUS_day7_plotfil$abs_eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_abs <- sort(quantile(AUS_day7_plotfil$abs_eff_size[which(AUS_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))
lenx <- length(s_nonsy_abs)
leny <- length(s_sy_abs)
if (leny < lenx)s_nonsy <- approx(1L:lenx, s_nonsy_abs, n = leny)$y
if (leny > lenx)s_sy <- approx(1L:leny, s_sy_abs, n = lenx)$y
p_AUS_sy_qq_abs<- ggplot() + 
  geom_point(aes(x=s_nonsy_abs, y=s_sy_abs)) +
  geom_abline(intercept = 0, slope = 1 ) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  style.print()+
  labs(title = 'AUS synon-nonsynon abs')

ggarrange(p_AUS_sy_qq,p_AUS_sy_qq_abs,ncol = 2,nrow = 1)


wilcox.test(AUS_day7_plotfil$abs_eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],AUS_day7_plotfil$abs_eff_size[which(AUS_day7_plotfil$group=='synonymous')],alternative ="greater")
wilcox.test(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group == 'synonymous')])


quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous')])
quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')])
wilcox.test(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous' & AUS_day7_plotfil$eff_size >quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],0.9))],
            AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous'& AUS_day7_plotfil$eff_size >quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous')],0.9))],alternative ="greater")

wilcox.test(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous' & AUS_day7_plotfil$eff_size <quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='nonsynonymous')],0.1))],
            AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous'& AUS_day7_plotfil$eff_size <quantile(AUS_day7_plotfil$eff_size[which(AUS_day7_plotfil$group=='synonymous')],0.1))],alternative ="less")

#TUU-----
TUU_day7_plotfil <- TUU_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_TUU-1)%>%
  mutate(abs_eff_size = abs(eff_size))%>%
  filter(eff_size != -1)

#TUU qq -----
s_nonsy_TUU <- sort(quantile(TUU_day7_plotfil$eff_size[which(TUU_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_TUU <- sort(quantile(TUU_day7_plotfil$eff_size[which(TUU_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_TUU_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_TUU, y=s_sy_TUU),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  "~ italic(P)< ~10^-12", x = -0.005, y = 0.017,size=50)


wilcox.test(TUU_sy_box$eff_size[TUU_sy_box$group=='synonymous'],TUU_sy_box$eff_size[TUU_sy_box$group=='nonsynonymous'],alternative = 'greater')
9.82e-13

#AUU-----
AUU_day7_plotfil <- AUU_day7_RSA %>%
  mutate(group= ifelse(AA_bef == AA_aft,'synonymous','nonsynonymous')) %>%
  mutate(type= ifelse(RSA > 0.2,'surfaced','buried'))%>%
  mutate(eff_size = fitness_AUU-1)%>%
  mutate(abs_eff_size = abs(eff_size))%>%
  filter(eff_size != -1)

#累积频率分布直方图（同义非同义）------
my_color = brewer.pal(n = 8, "Set2")[1:2]

figs2c_AUU <-ggplot(AUU_day7_plotfil, aes(x=eff_size, colour=group)) +
  scale_fill_manual(values= my_color)+
  geom_line(aes(y=..y..), stat="ecdf", size = 1)+
  theme_bw() + xlab('eff_size') +ylab('Cumulative frequency') +
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title = element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 0.5),
        legend.position = c(0.9,0.1),
        legend.justification  = c(0.9,0.1))+
  guides(colour=guide_legend(title=NULL))+
  labs(title = 'AUU synon-nonsynon')
#annotate('text', label = 'Wilcox test(abs) P = 0.0027(synon)', x = 0.03, y = 0.75,size=6)


# AUU qq------
s_nonsy_AUU <- sort(quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous')],probs = seq(0,1,0.01)))
s_sy_AUU <- sort(quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='synonymous')],probs = seq(0,1,0.01)))

p_AUU_sy_qq<- ggplot() + 
  geom_abline(intercept = 0, slope = 1,color = 'grey') +
  geom_hline(aes(yintercept=0), linetype = "dashed",color = 'grey') +
  geom_vline(aes(xintercept=0), linetype = "dashed",color = 'grey') +
  geom_point(aes(x=s_nonsy_AUU, y=s_sy_AUU),size=5) +
  style.print()+
  scale_x_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02),labels = c(-0.02,0,0.02))+
  scale_y_continuous(limits = c(-0.02,0.02),breaks = c(-0.02,0,0.02))+ 
  xlab('Nonsynonymous') +ylab('Synonymous') +
  theme(axis.text.x = element_text(size=100,hjust = 0.8),
        axis.title = element_text(size=120),axis.text.y=element_blank())+
  annotate('text',parse = TRUE,label =  "~ italic(P)< ~10^-60", x = -0.005, y = 0.017,size=50)


wilcox.test(AUU_sy_box$eff_size[AUU_sy_box$group=='synonymous'],AUU_sy_box$eff_size[AUU_sy_box$group=='nonsynonymous'],alternative = 'greater')$p.value
6.94219e-62

wilcox.test(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous')],AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group == 'synonymous')] )


quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='synonymous')])
quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous')])
wilcox.test(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous' & AUU_day7_plotfil$eff_size >quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous')],0.9))],
            AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='synonymous'& AUU_day7_plotfil$eff_size >quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='synonymous')],0.9))],alternative ="greater")

wilcox.test(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous' & AUU_day7_plotfil$eff_size <quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='nonsynonymous')],0.05))],
            AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='synonymous'& AUU_day7_plotfil$eff_size <quantile(AUU_day7_plotfil$eff_size[which(AUU_day7_plotfil$group=='synonymous')],0.05))],alternative ="less")


ggarrange(figs2c_AUS,figs2c_AUU,ncol = 2,nrow = 1)



#Figs2C sy nonsy 十组的box图--------

TGY_sy_box <- TGY_day7_plotfil[,c(12,14)] %>%
  mutate(type = 'TGY')
AGY_sy_box <- AGY_day7_plotfil[,c(12,14)] %>%
  mutate(type = 'AGY')
TGS_sy_box <- TGS_day7_plotfil[,c(12,14)] %>%
  mutate(type = 'TGS')
AGS_sy_box <- AGS_day7_plotfil[,c(12,14)] %>%
  mutate(type = 'AGS')
TUY_sy_box <- TUY_day7_plotfil[,c(4,14)] %>%
  mutate(type = 'TUY')
AUY_sy_box <- AUY_day7_plotfil[,c(4,14)] %>%
  mutate(type = 'AUY')
TUS_sy_box <- TUS_day7_plotfil[,c(4,14)] %>%
  mutate(type = 'TUS')
AUS_sy_box <- AUS_day7_plotfil[,c(4,14)] %>%
  mutate(type = 'AUS')
TUU_sy_box <- TUU_day7_plotfil[,c(4,14)] %>%
  mutate(type = 'TUU')
AUU_sy_box <- AUU_day7_plotfil[,c(4,14)] %>%
  mutate(type = 'AUU')

wilcox.test(AGY_sy_box$eff_size[AGY_sy_box$group=='synonymous'],AGY_sy_box$eff_size[AGY_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.2147
wilcox.test(TGY_sy_box$eff_size[TGY_sy_box$group=='synonymous'],TGY_sy_box$eff_size[TGY_sy_box$group=='nonsynonymous'],alternative = 'greater')
2.2e-16
wilcox.test(TGS_sy_box$eff_size[TGS_sy_box$group=='synonymous'],TGS_sy_box$eff_size[TGS_sy_box$group=='nonsynonymous'],alternative = 'greater')
4.885e-15
wilcox.test(AGS_sy_box$eff_size[AGS_sy_box$group=='synonymous'],AGS_sy_box$eff_size[AGS_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.06373
wilcox.test(AUY_sy_box$eff_size[AUY_sy_box$group=='synonymous'],AUY_sy_box$eff_size[AUY_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.01837
wilcox.test(TUY_sy_box$eff_size[TUY_sy_box$group=='synonymous'],TUY_sy_box$eff_size[TUY_sy_box$group=='nonsynonymous'],alternative = 'greater')
0.01001
wilcox.test(TUS_sy_box$eff_size[TUS_sy_box$group=='synonymous'],TUS_sy_box$eff_size[TUS_sy_box$group=='nonsynonymous'],alternative = 'greater')
2.2e-16
wilcox.test(AUS_sy_box$eff_size[AUS_sy_box$group=='synonymous'],AUS_sy_box$eff_size[AUS_sy_box$group=='nonsynonymous'],alternative = 'greater')
2.383e-12
wilcox.test(TUU_sy_box$eff_size[TUU_sy_box$group=='synonymous'],TUU_sy_box$eff_size[TUU_sy_box$group=='nonsynonymous'],alternative = 'greater')
9.82e-13
wilcox.test(AUU_sy_box$eff_size[AUU_sy_box$group=='synonymous'],AUU_sy_box$eff_size[AUU_sy_box$group=='nonsynonymous'],alternative = 'greater')
2.2e-16

TA_UY_paired_data <-merge(TUY_sy_box,AUY_sy_box,by='mut',all=T)
TA_UY_paired_data <- na.omit(TA_UY_paired_data)
wilcox.test(TA_UY_paired_data$eff_size.x,TA_UY_paired_data$eff_size.y,paired = TRUE,alternative = "greater")
median(TA_UY_paired_data$eff_size.x)
median(TA_UY_paired_data$eff_size.y)
mean(TA_UY_paired_data$eff_size.x)
mean(TA_UY_paired_data$eff_size.y)

Figs2c_10groupbox <- rbind(TGY_sy_box,AGY_sy_box,TGS_sy_box,AGS_sy_box,TUY_sy_box,AUY_sy_box,TUS_sy_box,AUS_sy_box,TUU_sy_box,AUU_sy_box)


library(RColorBrewer)
my_orange = brewer.pal(n = 8, "Set2")[3:6]

Figs2c_10groupbox %>%
  mutate(x = factor(type,levels = c('TGY','AGY','TGS','AGS','TUY','AUY','TUS','AUS','TUU','AUU'),order=T))%>%
  mutate(set = factor(group,levels = c('synonymous','nonsynonymous'),order=T))%>%
  ggplot(aes(x = x, y = eff_size,fill=set)) +
  scale_fill_manual(values = my_orange)+
  ylab('eff_size') + xlab('') +
  geom_boxplot(alpha=0.8,width=0.6,scale = 'width')+
  #facet_grid(.~set, scales = "free", switch = "x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ ylim(-0.02,0.02) +
  style.print()


# FigS3 c)换成了qq图合集--------------

pdf("~/figS3c_GFP.pdf",width = 17,height = 15)
ggarrange(p_TGY_sy_qq,p_TGS_sy_qq,p_AGY_sy_qq,p_AGS_sy_qq)
dev.off()

pdf("~/figS3c_URA3.pdf",width = 23,height = 14)
ggarrange(p_TUY_sy_qq,p_TUS_sy_qq,p_TUU_sy_qq,p_AUY_sy_qq,p_AUS_sy_qq,p_AUU_sy_qq)
dev.off()


pdf("~/figS3c_GFP_URA3.pdf",width = 79,height = 29)
ggarrange(p_TGY_sy_qq,p_TGS_sy_qq,p_TUY_sy_qq,p_TUS_sy_qq,p_TUU_sy_qq,p_AGY_sy_qq,p_AGS_sy_qq,p_AUY_sy_qq,p_AUS_sy_qq,p_AUU_sy_qq,ncol = 5,nrow = 2)
dev.off()

pdf("~/figS3c_GFP_URA3_noGS.pdf",width = 63,height = 29)
ggarrange(p_TGY_sy_qq,p_TUY_sy_qq,p_TUS_sy_qq,p_TUU_sy_qq,p_AGY_sy_qq,p_AUY_sy_qq,p_AUS_sy_qq,p_AUU_sy_qq,ncol = 4,nrow = 2)
dev.off()


# FigS4 AB------
#GFP-----

#外部misfo 0.2
TGY_day7_surfaced_misfolding <- filter(TGY_day7_misfolding, RSA>0.2)

p_TGY_misfo_surfaced<- ggplot(TGY_day7_surfaced_misfolding,aes(TGY_day7_surfaced_misfolding$misfolding,TGY_day7_surfaced_misfolding$fitness_TGY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.09  P=0.006', x = 57, y = 1.02, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
    panel.grid =element_blank(),
  panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TGY surface RSA0.2 misfolding")

cor.test(TGY_day7_surfaced_misfolding$misfolding,TGY_day7_surfaced_misfolding$fitness_TGY,method = "s")


AGY_day7_surfaced_misfolding <- filter(AGY_day7_misfolding, RSA>0.2)
p_AGY_misfo_surfaced<- ggplot(AGY_day7_surfaced_misfolding,aes(AGY_day7_surfaced_misfolding$misfolding,AGY_day7_surfaced_misfolding$fitness_AGY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.03  P=0.326', x = 56, y = 1.02, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AGY surface RSA0.2 misfolding")

cor.test(AGY_day7_surfaced_misfolding$misfolding,AGY_day7_surfaced_misfolding$fitness_AGY,method = "s")


TGS_day7_surfaced_misfolding <- filter(TGS_day7_misfolding, RSA>0.2)
p_TGS_misfo_surfaced<- ggplot(TGS_day7_surfaced_misfolding,aes(TGS_day7_surfaced_misfolding$misfolding,TGS_day7_surfaced_misfolding$fitness_TGS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.06  P=0.051', x = 56, y = 1.06, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TGS surface RSA0.2 misfolding")

cor.test(TGS_day7_surfaced_misfolding$misfolding,TGS_day7_surfaced_misfolding$fitness_TGS,method = "s")


AGS_day7_surfaced_misfolding <- filter(AGS_day7_misfolding, RSA>0.2)
p_AGS_misfo_surfaced<- ggplot(AGS_day7_surfaced_misfolding,aes(AGS_day7_surfaced_misfolding$misfolding,AGS_day7_surfaced_misfolding$fitness_AGS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.004  P=0.913', x = 22, y = 1.08, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ ylim(0.96,1.08) + xlim(0,30)+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AGS surface RSA0.2 misfolding")

cor.test(AGS_day7_surfaced_misfolding$misfolding,AGS_day7_surfaced_misfolding$fitness_AGS,method = "s")

ggarrange(p_TGY_misfo_surfaced,p_TGS_misfo_surfaced,p_AGY_misfo_surfaced,p_AGS_misfo_surfaced,ncol = 2,nrow = 2)


#外部misfo 0.4
TGY_day7_surfaced_misfolding <- filter(TGY_day7_misfolding, RSA>0.5)

p_TGY_misfo_surfaced<- ggplot(TGY_day7_surfaced_misfolding,aes(TGY_day7_surfaced_misfolding$misfolding,TGY_day7_surfaced_misfolding$fitness_TGY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.074  P=0.225', x = 9, y = 1.01, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TGY surface RSA0.5 misfolding")

cor.test(TGY_day7_surfaced_misfolding$misfolding,TGY_day7_surfaced_misfolding$fitness_TGY,method = "s")


AGY_day7_surfaced_misfolding <- filter(AGY_day7_misfolding, RSA>0.5)
p_AGY_misfo_surfaced<- ggplot(AGY_day7_surfaced_misfolding,aes(AGY_day7_surfaced_misfolding$misfolding,AGY_day7_surfaced_misfolding$fitness_AGY)) +
  geom_point() +
  annotate('text', label = 'ρ=0.02  P=0.78', x = 9, y = 1.02, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AGY surface RSA0.5 misfolding")

cor.test(AGY_day7_surfaced_misfolding$misfolding,AGY_day7_surfaced_misfolding$fitness_AGY,method = "s")


TGS_day7_surfaced_misfolding <- filter(TGS_day7_misfolding, RSA>0.5)
p_TGS_misfo_surfaced<- ggplot(TGS_day7_surfaced_misfolding,aes(TGS_day7_surfaced_misfolding$misfolding,TGS_day7_surfaced_misfolding$fitness_TGS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.06  P=0.336', x = 9, y = 1.06, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TGS surface RSA0.5 misfolding")

cor.test(TGS_day7_surfaced_misfolding$misfolding,TGS_day7_surfaced_misfolding$fitness_TGS,method = "s")


AGS_day7_surfaced_misfolding <- filter(AGS_day7_misfolding, RSA>0.5)
p_AGS_misfo_surfaced<- ggplot(AGS_day7_surfaced_misfolding,aes(AGS_day7_surfaced_misfolding$misfolding,AGS_day7_surfaced_misfolding$fitness_AGS)) +
  geom_point() +
  annotate('text', label = 'ρ=0.04  P=0.476', x = 9, y = 1.08, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ ylim(0.96,1.08)+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AGS surface RSA0.5 misfolding")

cor.test(AGS_day7_surfaced_misfolding$misfolding,AGS_day7_surfaced_misfolding$fitness_AGS,method = "s")

ggarrange(p_TGY_misfo_surfaced,p_TGS_misfo_surfaced,p_AGY_misfo_surfaced,p_AGS_misfo_surfaced,ncol = 2,nrow = 2)





#内部misinteraction
TGY_day7_buried <- filter(TGY_day7_geno, RSA < 0.2)
TGY_day7_buried <- TGY_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TGY=mean(fitness_TGY))
cor.test(TGY_day7_buried$hyb_t,TGY_day7_buried$fitness_TGY,method = "s")

TGY_day7_buried$hyb_RSA <- TGY_day7_buried$hyb_t * TGY_day7_buried$RSA
p_TGY_misin_buried <-ggplot(TGY_day7_buried,aes(TGY_day7_buried$hyb_RSA,TGY_day7_buried$fitness_TGY)) +
  annotate('text', label = 'ρ=-0.03  P=0.391', x = 1, y = 1.015, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) +
  labs(title = "TGY buried misinteraction")

cor.test(TGY_day7_buried$hyb_RSA,TGY_day7_buried$fitness_TGY,method = "s")

TGS_day7_buried <- filter(TGS_day7_geno, RSA < 0.2)
TGS_day7_buried <- TGS_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TGS=mean(fitness_TGS))
cor.test(TGS_day7_buried$hyb_t,TGS_day7_buried$fitness_TGS,method = "s")

TGS_day7_buried$hyb_RSA <- TGS_day7_buried$hyb_t * TGS_day7_buried$RSA
p_TGS_misin_buried <-ggplot(TGS_day7_buried,aes(TGS_day7_buried$hyb_RSA,TGS_day7_buried$fitness_TGS)) +
  annotate('text', label = 'ρ=-0.004  P=0.91', x = 1, y = 1.057, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TGS buried misinteraction")
cor.test(TGS_day7_buried$hyb_RSA,TGS_day7_buried$fitness_TGS,method = "s")


AGY_day7_buried <- filter(AGY_day7_geno, RSA < 0.2)
AGY_day7_buried <- AGY_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AGY=mean(fitness_AGY))
cor.test(AGY_day7_buried$hyb_t,AGY_day7_buried$fitness_AGY,method = "s")

AGY_day7_buried$hyb_RSA <- AGY_day7_buried$hyb_t * AGY_day7_buried$RSA
p_AGY_misin_buried <-ggplot(AGY_day7_buried,aes(AGY_day7_buried$hyb_RSA,AGY_day7_buried$fitness_AGY)) +
  annotate('text', label = 'ρ=-0.09  P=0.009', x = 1, y = 1.023, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AGY buried misinteraction")
cor.test(AGY_day7_buried$hyb_RSA,AGY_day7_buried$fitness_AGY,method = "s")



AGS_day7_buried <- filter(AGS_day7_geno, RSA < 0.2)
AGS_day7_buried <- AGS_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AGS=mean(fitness_AGS))
cor.test(AGS_day7_buried$hyb_t,AGS_day7_buried$fitness_AGS,method = "s")

AGS_day7_buried$hyb_RSA <- AGS_day7_buried$hyb_t * AGS_day7_buried$RSA
p_AGS_misin_buried <-ggplot(AGS_day7_buried,aes(AGS_day7_buried$hyb_RSA,AGS_day7_buried$fitness_AGS)) +
  annotate('text', label = 'ρ=-0.017  P=0.58', x = 1, y = 1.075, size=7, colour ='#000000')+
  geom_point() + ylim(0.96,1.08)+
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AGS buried misinteraction")
cor.test(AGS_day7_buried$hyb_RSA,AGS_day7_buried$fitness_AGS,method = "s")

ggarrange(p_TGY_misin_buried,p_TGS_misin_buried,p_AGY_misin_buried,p_AGS_misin_buried,ncol = 2,nrow = 2)



#URA3-----

#RSA0.2 外部misfolding-----
TUY_day7_surfaced_misfolding <- filter(TUY_day7_misfolding, RSA>0.2)

cor.test(TUY_day7_surfaced_misfolding$misfolding,TUY_day7_surfaced_misfolding$fitness_TUY,method = "s")
p_TUY_misfo_surfaced<- ggplot(TUY_day7_surfaced_misfolding,aes(TUY_day7_surfaced_misfolding$misfolding,TUY_day7_surfaced_misfolding$fitness_TUY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.06 P=0.049', x = 40, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUY surface RSA0.2 misfolding")


AUY_day7_surfaced_misfolding <- filter(AUY_day7_misfolding, RSA>0.2)
cor.test(AUY_day7_surfaced_misfolding$misfolding,AUY_day7_surfaced_misfolding$fitness_AUY,method = "s")
p_AUY_misfo_surfaced<- ggplot(AUY_day7_surfaced_misfolding,aes(AUY_day7_surfaced_misfolding$misfolding,AUY_day7_surfaced_misfolding$fitness_AUY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.05  P=0.1125', x = 40, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUY surface RSA0.2 misfolding")



TUS_day7_surfaced_misfolding <- filter(TUS_day7_misfolding, RSA>0.2)
cor.test(TUS_day7_surfaced_misfolding$misfolding,TUS_day7_surfaced_misfolding$fitness_TUS,method = "s")
p_TUS_misfo_surfaced<- ggplot(TUS_day7_surfaced_misfolding,aes(TUS_day7_surfaced_misfolding$misfolding,TUS_day7_surfaced_misfolding$fitness_TUS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.113  P=0.0003', x = 30, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUS surface RSA0.2 misfolding")



AUS_day7_surfaced_misfolding <- filter(AUS_day7_misfolding, RSA>0.2)
cor.test(AUS_day7_surfaced_misfolding$misfolding,AUS_day7_surfaced_misfolding$fitness_AUS,method = "s")
p_AUS_misfo_surfaced<- ggplot(AUS_day7_surfaced_misfolding,aes(AUS_day7_surfaced_misfolding$misfolding,AUS_day7_surfaced_misfolding$fitness_AUS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.085  P=0.006', x = 20, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ xlim(0,30)+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUS surface RSA0.2 misfolding")


ggarrange(p_TUY_misfo_surfaced,p_TUS_misfo_surfaced,p_AUY_misfo_surfaced,p_AUS_misfo_surfaced,ncol = 2,nrow = 2)


TUU_day7_surfaced_misfolding <- filter(TUU_day7_misfolding, RSA>0.2)
cor.test(TUU_day7_surfaced_misfolding$misfolding,TUU_day7_surfaced_misfolding$fitness_TUU,method = "s")
p_TUU_misfo_surfaced<- ggplot(TUU_day7_surfaced_misfolding,aes(TUU_day7_surfaced_misfolding$misfolding,TUU_day7_surfaced_misfolding$fitness_TUU)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.113  P=0.0003', x = 30, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUU surface RSA0.2 misfolding")



AUU_day7_surfaced_misfolding <- filter(AUU_day7_misfolding, RSA>0.2)
cor.test(AUU_day7_surfaced_misfolding$misfolding,AUU_day7_surfaced_misfolding$fitness_AUU,method = "s")
p_AUU_misfo_surfaced<- ggplot(AUU_day7_surfaced_misfolding,aes(AUU_day7_surfaced_misfolding$misfolding,AUU_day7_surfaced_misfolding$fitness_AUU)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.085  P=0.006', x = 20, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ xlim(0,30)+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUU surface RSA0.2 misfolding")






#RSA0.4 外部misfolding-------
TUY_day7_surfaced_misfolding <- filter(TUY_day7_misfolding, RSA>0.4)

p_TUY_misfo_surfaced<- ggplot(TUY_day7_surfaced_misfolding,aes(TUY_day7_surfaced_misfolding$misfolding,TUY_day7_surfaced_misfolding$fitness_TUY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.031 P=0.5', x = 20, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUY surface RSA0.4 misfolding")

cor.test(TUY_day7_surfaced_misfolding$misfolding,TUY_day7_surfaced_misfolding$fitness_TUY,method = "s")


AUY_day7_surfaced_misfolding <- filter(AUY_day7_misfolding, RSA>0.4)
p_AUY_misfo_surfaced<- ggplot(AUY_day7_surfaced_misfolding,aes(AUY_day7_surfaced_misfolding$misfolding,AUY_day7_surfaced_misfolding$fitness_AUY)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.032  P=0.483', x = 20, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
   labs(title = "AUY surface RSA0.4 misfolding")
cor.test(AUY_day7_surfaced_misfolding$misfolding,AUY_day7_surfaced_misfolding$fitness_AUY,method = "s")


TUS_day7_surfaced_misfolding <- filter(TUS_day7_misfolding, RSA>0.4)
p_TUS_misfo_surfaced<- ggplot(TUS_day7_surfaced_misfolding,aes(TUS_day7_surfaced_misfolding$misfolding,TUS_day7_surfaced_misfolding$fitness_TUS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.069  P=0.126', x = 20, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUS surface RSA0.4 misfolding")

cor.test(TUS_day7_surfaced_misfolding$misfolding,TUS_day7_surfaced_misfolding$fitness_TUS,method = "s")


AUS_day7_surfaced_misfolding <- filter(AUS_day7_misfolding, RSA>0.4)
p_AUS_misfo_surfaced<- ggplot(AUS_day7_surfaced_misfolding,aes(AUS_day7_surfaced_misfolding$misfolding,AUS_day7_surfaced_misfolding$fitness_AUS)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.094  P=0.038', x = 22, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ ylim(0.96,1.08) + xlim(0,30)+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUS surface RSA0.4 misfolding")

cor.test(AUS_day7_surfaced_misfolding$misfolding,AUS_day7_surfaced_misfolding$fitness_AUS,method = "s")

ggarrange(p_TUY_misfo_surfaced,p_TUS_misfo_surfaced,p_AUY_misfo_surfaced,p_AUS_misfo_surfaced,ncol = 2,nrow = 2)

TUU_day7_surfaced_misfolding <- filter(TUU_day7_misfolding, RSA>0.4)
cor.test(TUU_day7_surfaced_misfolding$misfolding,TUU_day7_surfaced_misfolding$fitness_TUU,method = "s")
p_TUU_misfo_surfaced<- ggplot(TUU_day7_surfaced_misfolding,aes(TUU_day7_surfaced_misfolding$misfolding,TUU_day7_surfaced_misfolding$fitness_TUU)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.113  P=0.0003', x = 30, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ 
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUU surface RSA0.4 misfolding")



AUU_day7_surfaced_misfolding <- filter(AUU_day7_misfolding, RSA>0.4)
cor.test(AUU_day7_surfaced_misfolding$misfolding,AUU_day7_surfaced_misfolding$fitness_AUU,method = "s")
p_AUU_misfo_surfaced<- ggplot(AUU_day7_surfaced_misfolding,aes(AUU_day7_surfaced_misfolding$misfolding,AUU_day7_surfaced_misfolding$fitness_AUU)) +
  geom_point() +
  annotate('text', label = 'ρ=-0.085  P=0.006', x = 20, y = 1.05, size=7, colour ='#000000')+
  geom_smooth(method = 'lm',se=F) + xlab("P_misfolding") + ylab("fitness")+ xlim(0,30)+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUU surface RSA0.4 misfolding")

#内部misinteraction
TUY_day7_buried <- filter(TUY_day7_geno, RSA < 0.2)
TUY_day7_buried <- TUY_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TUY=mean(fitness_TUY))
cor.test(TUY_day7_buried$hyb_t,TUY_day7_buried$fitness_TUY,method = "s")

TUY_day7_buried$hyb_RSA <- TUY_day7_buried$hyb_t * TUY_day7_buried$RSA
p_TUY_misin_buried <-ggplot(TUY_day7_buried,aes(TUY_day7_buried$hyb_RSA,TUY_day7_buried$fitness_TUY)) +
  annotate('text', label = 'ρ=0.055  P=0.106', x = 0.5, y = 1.05, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) +
  labs(title = "TUY buried misinteraction")

cor.test(TUY_day7_buried$hyb_RSA,TUY_day7_buried$fitness_TUY,method = "s")


TUS_day7_buried <- filter(TUS_day7_geno, RSA < 0.2)
TUS_day7_buried <- TUS_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TUS=mean(fitness_TUS))
cor.test(TUS_day7_buried$hyb_t,TUS_day7_buried$fitness_TUS,method = "s")

TUS_day7_buried$hyb_RSA <- TUS_day7_buried$hyb_t * TUS_day7_buried$RSA
p_TUS_misin_buried <-ggplot(TUS_day7_buried,aes(TUS_day7_buried$hyb_RSA,TUS_day7_buried$fitness_TUS)) +
  annotate('text', label = 'ρ=0.0065  P=0.846', x = 0.2, y = 1.05, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "TUS buried misinteraction")
cor.test(TUS_day7_buried$hyb_RSA,TUS_day7_buried$fitness_TUS,method = "s")


AUY_day7_buried <- filter(AUY_day7_geno, RSA < 0.2)
AUY_day7_buried <- AUY_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AUY=mean(fitness_AUY))
cor.test(AUY_day7_buried$hyb_t,AUY_day7_buried$fitness_AUY,method = "s")

AUY_day7_buried$hyb_RSA <- AUY_day7_buried$hyb_t * AUY_day7_buried$RSA
p_AUY_misin_buried <-ggplot(AUY_day7_buried,aes(AUY_day7_buried$hyb_RSA,AUY_day7_buried$fitness_AUY)) +
  annotate('text', label = 'ρ=0.07  P=0.029', x = 0.5, y = 1.05, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUY buried misinteraction")
cor.test(AUY_day7_buried$hyb_RSA,AUY_day7_buried$fitness_AUY,method = "s")



AUS_day7_buried <- filter(AUS_day7_geno, RSA < 0.2)
AUS_day7_buried <- AUS_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AUS=mean(fitness_AUS))
cor.test(AUS_day7_buried$hyb_t,AUS_day7_buried$fitness_AUS,method = "s")

AUS_day7_buried$hyb_RSA <- AUS_day7_buried$hyb_t * AUS_day7_buried$RSA
p_AUS_misin_buried <-ggplot(AUS_day7_buried,aes(AUS_day7_buried$hyb_RSA,AUS_day7_buried$fitness_AUS)) +
  annotate('text', label = 'ρ=-0.034  P=0.31', x = 0.5, y = 1.06, size=7, colour ='#000000')+
  geom_point() + ylim(0.96,1.08)+
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUS buried misinteraction")
cor.test(AUS_day7_buried$hyb_RSA,AUS_day7_buried$fitness_AUS,method = "s")

ggarrange(p_TUY_misin_buried,p_TUS_misin_buried,p_AUY_misin_buried,p_AUS_misin_buried,ncol = 2,nrow = 2)


TUU_day7_buried <- filter(TUU_day7_geno, RSA < 0.2)
TUU_day7_buried <- TUU_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_TUU=mean(fitness_TUU))
cor.test(TUU_day7_buried$hyb_t,TUU_day7_buried$fitness_TUU,method = "s")

TUU_day7_buried$hyb_RSA <- TUU_day7_buried$hyb_t * TUU_day7_buried$RSA
p_TUU_misin_buried <-ggplot(TUU_day7_buried,aes(TUU_day7_buried$hyb_RSA,TUU_day7_buried$fitness_TUU)) +
  annotate('text', label = 'ρ=0.055  P=0.106', x = 0.5, y = 1.05, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) +
  labs(title = "TUU buried misinteraction")

cor.test(TUU_day7_buried$hyb_RSA,TUU_day7_buried$fitness_TUU,method = "s")


AUU_day7_buried <- filter(AUU_day7_geno, RSA < 0.2)
AUU_day7_buried <- AUU_day7_buried %>%
  dplyr::group_by(AA_pos,RSA,hyb_t) %>%
  dplyr::summarise(fitness_AUU=mean(fitness_AUU))
cor.test(AUU_day7_buried$hyb_t,AUU_day7_buried$fitness_AUU,method = "s")

AUU_day7_buried$hyb_RSA <- AUU_day7_buried$hyb_t * AUU_day7_buried$RSA
p_AUU_misin_buried <-ggplot(AUU_day7_buried,aes(AUU_day7_buried$hyb_RSA,AUU_day7_buried$fitness_AUU)) +
  annotate('text', label = 'ρ=0.07  P=0.029', x = 0.5, y = 1.05, size=7, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("Diff_hydrophobicity*RSA") + ylab("fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=15),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5))+
  labs(title = "AUU buried misinteraction")
cor.test(AUU_day7_buried$hyb_RSA,AUU_day7_buried$fitness_AUU,method = "s")



# misinteraction 全局所有样本 ----------

#TGY-------

TGY_day7_formisin <- TGY_day7_geno

TGY_day7_formisin$hyb_RSA <- TGY_day7_formisin$hyb_t * TGY_day7_formisin$RSA
p_TGY_misin_all <- ggplot(TGY_day7_formisin,aes(TGY_day7_formisin$hyb_RSA,TGY_day7_formisin$fitness_TGY)) +
  annotate('text', label =  "'ρ= -0.08 ' ~ italic(P)< ~10^-3", parse = T, x = 3.3, y = 1.0155, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(TGY_day7_formisin$hyb_RSA,TGY_day7_formisin$fitness_TGY,method = "s")


#AGY-------

AGY_day7_formisin <- AGY_day7_geno

AGY_day7_formisin$hyb_RSA <- AGY_day7_formisin$hyb_t * AGY_day7_formisin$RSA
p_AGY_misin_all <- ggplot(AGY_day7_formisin,aes(AGY_day7_formisin$hyb_RSA,AGY_day7_formisin$fitness_AGY)) +
  annotate('text', label =  "'ρ= -0.05 ' ~ italic(P)< ~0.05", parse = T, x = 3.3, y = 1.027, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(AGY_day7_formisin$hyb_RSA,AGY_day7_formisin$fitness_AGY,method = "s")


#TUY-------

TUY_day7_formisin <- TUY_day7_geno

TUY_day7_formisin$hyb_RSA <- TUY_day7_formisin$hyb_t * TUY_day7_formisin$RSA
TUY_day7_formisin <- TUY_day7_formisin[-which(TUY_day7_formisin$fitness_TUY == 0),]
p_TUY_misin_all <- ggplot(TUY_day7_formisin,aes(TUY_day7_formisin$hyb_RSA,TUY_day7_formisin$fitness_TUY)) +
  annotate('text', label =  "'ρ= 0.04 ' ~ italic(P) == 0.06", parse = T, x = 3.7, y = 1.053, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(TUY_day7_formisin$hyb_RSA,TUY_day7_formisin$fitness_TUY,method = "s")


#AUY-------

AUY_day7_formisin <- AUY_day7_geno

AUY_day7_formisin$hyb_RSA <- AUY_day7_formisin$hyb_t * AUY_day7_formisin$RSA
AUY_day7_formisin <- AUY_day7_formisin[-which(AUY_day7_formisin$fitness_AUY == 0),]
p_AUY_misin_all <- ggplot(AUY_day7_formisin,aes(AUY_day7_formisin$hyb_RSA,AUY_day7_formisin$fitness_AUY)) +
  annotate('text', label =  "'ρ= 0.01 ' ~ italic(P) == 0.77", parse = T, x = 3.7, y = 1.057, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(AUY_day7_formisin$hyb_RSA,AUY_day7_formisin$fitness_AUY,method = "s")



#TUS-------

TUS_day7_formisin <- TUS_day7_geno

TUS_day7_formisin$hyb_RSA <- TUS_day7_formisin$hyb_t * TUS_day7_formisin$RSA
TUS_day7_formisin <- TUS_day7_formisin[-which(TUS_day7_formisin$fitness_TUS == 0),]
p_TUS_misin_all <- ggplot(TUS_day7_formisin,aes(TUS_day7_formisin$hyb_RSA,TUS_day7_formisin$fitness_TUS)) +
  annotate('text', label =  "'ρ= 0.03 ' ~ italic(P) == 0.08", parse = T, x = 3.7, y = 1.052, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(TUS_day7_formisin$hyb_RSA,TUS_day7_formisin$fitness_TUS,method = "s")

#AUS-------

AUS_day7_formisin <- AUS_day7_geno

AUS_day7_formisin$hyb_RSA <- AUS_day7_formisin$hyb_t * AUS_day7_formisin$RSA
AUS_day7_formisin <- AUS_day7_formisin[-which(AUS_day7_formisin$fitness_AUS == 0),]
p_AUS_misin_all <- ggplot(AUS_day7_formisin,aes(AUS_day7_formisin$hyb_RSA,AUS_day7_formisin$fitness_AUS)) +
  annotate('text', label =  "'ρ= 0.02 ' ~ italic(P) == 0.35", parse = T, x = 3.7, y = 1.056, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(AUS_day7_formisin$hyb_RSA,AUS_day7_formisin$fitness_AUS,method = "s")





#TUU-------

TUU_day7_formisin <- TUU_day7_geno

TUU_day7_formisin$hyb_RSA <- TUU_day7_formisin$hyb_t * TUU_day7_formisin$RSA
TUU_day7_formisin <- TUU_day7_formisin[-which(TUU_day7_formisin$fitness_TUU == 0),]
p_TUU_misin_all <- ggplot(TUU_day7_formisin,aes(TUU_day7_formisin$hyb_RSA,TUU_day7_formisin$fitness_TUU)) +
  annotate('text', label =  "'ρ= 0.003 ' ~ italic(P) == 0.9", parse = T, x = 3.7, y = 1.055, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(TUU_day7_formisin$hyb_RSA,TUU_day7_formisin$fitness_TUU,method = "s")



#AUY-------

AUU_day7_formisin <- AUU_day7_geno

AUU_day7_formisin$hyb_RSA <- AUU_day7_formisin$hyb_t * AUU_day7_formisin$RSA
AUU_day7_formisin <- AUU_day7_formisin[-which(AUU_day7_formisin$fitness_AUU == 0),]
p_AUU_misin_all <- ggplot(AUU_day7_formisin,aes(AUU_day7_formisin$hyb_RSA,AUU_day7_formisin$fitness_AUU)) +
  annotate('text', label =  "'ρ= 0.003 ' ~ italic(P) == 0.89", parse = T, x = 3.7, y = 1.046, size=25, colour ='#000000')+
  geom_point() +
  geom_smooth(method = 'lm',se=F) + xlab("RSA*Hydrophobicity") + ylab("Fitness")+
  theme_bw() + 
  theme(axis.text=element_text(size=70),
        axis.title.x =element_text(size=70),
        axis.title.y=element_text(size=70),
        panel.grid =element_blank(),
        panel.border= element_rect(color = 'black', size = 1.5)) 

cor.test(AUU_day7_formisin$hyb_RSA,AUU_day7_formisin$fitness_AUU,method = "s")



cairo_pdf("~/figS4c.pdf",width = 70,height = 29)
ggarrange(p_TGY_misin_all,p_TUY_misin_all,p_TUS_misin_all,p_TUU_misin_all,p_AGY_misin_all,p_AUY_misin_all,p_AUS_misin_all,p_AUU_misin_all,ncol = 4,nrow = 2)
dev.off()





