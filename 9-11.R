area_fitness<- read.csv("9-12_dis10_average_fitness.csv",header = F)
fitness_data<-as.data.frame(area_fitness)
names(fitness_data)<-c("pos","fitness")
fitness_data$pos <- fitness_data$pos -2


p <- ggplot(fitness_data,aes(pos,fitness)) +
  geom_point() +
  geom_smooth(method="lm") + xlab("pos") + ylab("effective size_10") 
p + scale_x_continuous(breaks=seq(0, 270, 10))


cutoff_npass <- read.csv("10-16_reads_npass.csv",header = T)
names(cutoff_npass) <- c("cutoff","type","value")

cutoff_npass$type <-factor(cutoff_npass$type, levels = c("umi","geno"))

options(scipen=200)

ggplot(cutoff_npass,aes(x=cutoff,y=value,fill=type))+ geom_bar(stat = "identity",position = "dodge")+ 
  scale_x_continuous(breaks=seq(0, 30, 5))+ scale_y_continuous(expand =c(0,0),breaks=seq(0,90000,10000))+ style.print() +
   xlab("cutoff > x") +ylab("count") + theme(legend.title = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                             axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) 


reads_cutoff <- read.csv("11-6cutoff_reads.csv",header = T)

options(scipen=1)

ggplot(reads_cutoff,aes(x=cutoff,y=count,fill=type))+ geom_bar(stat = "identity",position = "dodge",width = 0.6)+ 
  scale_x_continuous(breaks=seq(0, 10, 1))+ 
  scale_y_continuous(expand =c(0,0),breaks=seq(0, 3.5, 0.5)) + style.print() +
  xlab("Number of npass >= x") +ylab("Number of reads") + theme(legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                            axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))

# npass_cell number 分链了

npass_cellnum <- read.csv("11-6npass_num_density.csv",header = T)

npass_cellnum$npass <- factor(npass_cellnum$npass, levels=unique(npass_cellnum$npass))

ggplot(npass_cellnum,aes(x=npass,y=cell))+ geom_bar(stat = "identity",position = "dodge",width = 0.6,color = "grey")+ 
  scale_y_continuous(expand =c(0.1,0),breaks=seq(0, 500000, 25000)) + style.print() +
  xlab("npass = x") +ylab("Number of zwm") + theme(legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                       axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))

#没分链

npass_cellnum_nosplit <- read.csv("11-7detech_notsplit.csv",header = T)

npass_cellnum_nosplit$npass <- factor(npass_cellnum_nosplit$npass, levels=unique(npass_cellnum_nosplit$npass))

ggplot(npass_cellnum_nosplit,aes(x=npass,y=cell))+ geom_bar(stat = "identity",position = "dodge",width = 0.6,color = "grey")+ 
  scale_y_continuous(expand =c(0.1,0),breaks=seq(0, 500000, 25000)) + style.print() +
  xlab("npass = x") +ylab("Number of zwm") + theme(legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                   axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))



# ,labels = function(x) format(x, scientific = TRUE)

nofilter_reads <- read.csv("10-16_nofilter_reads.csv",header = T)
names(nofilter_reads) <- c("cutoff","type","value")
nofilter_reads$type <-factor(nofilter_reads$type, levels = c("geno","umi","reads"))
ggplot(nofilter_reads,aes(x=cutoff,y=value,fill=type))+ geom_bar(stat = "identity",position = "stack")+ 
  scale_x_continuous(breaks=seq(0, 30, 5))+ scale_y_continuous(breaks=seq(0, 700000, 50000)) + style.print()


#单点突变组合图

cell_combine <- read.csv("10-30_singlemut_combine.csv",header = T)
names(cell_combine) <- c("cell_num","type","value")

cell_combine$cell_num <-factor(cell_combine$cell_num, levels = c("single_cell","two_cells","three_cells"))


ggplot(cell_combine,aes(x=cell_num,y=value,fill=type))+ geom_bar(stat = "identity",position = "dodge")+ 
   scale_y_continuous(expand =c(0,0),breaks=seq(0,12000,1000))+ style.print() +
  xlab("Number of cells") +ylab("count") + theme(legend.title = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                            axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) 



#11-13 画质量图

npass_quality <- read.csv("54136f_npass_quality.csv",header = F)
names(npass_quality) <- c("npass","quality")

npass_quality$npass <- factor(npass_quality$npass, levels=unique(npass_quality$npass))

ggplot(npass_quality,aes(x=npass,y=quality))+ geom_bar(stat = "identity",position = "dodge",width = 0.4,color = "grey")+ 
  scale_y_continuous(expand =c(0.1,0),breaks=seq(0, 140, 20)) + style.print() +
  xlab("npass = x") +ylab("Average quality") + theme(legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                                                   axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 18))





