sub(".*(TTCAACCAAGTTG.*)","\\1",seq,perl=T);
source("https://bioconductor.org/biocLite.R")

sub(".*(TTCAACCAAG.*)","\\1",seq,perl=T)

source("/mnt/data/home/phil/lib/R/biocLite.R")
biocLite("seqinr")

dfDouble<-subset(dfDouble,select=-mutPos1)

me<-mclapply(1:10000,function(i){
  
  tt<-mean(as.numeric(te[i,]))
  new<-data.frame(mean=tt)
},mc.cores = 30)%>%rbind.fill()


apply(te,1,mean)
a<-sub("AGTGTATACAAATTTTAAAGTGACTCTTAGGTTTTAAAACGAAAATTCTTATTCTTGAGTAACTCTTTCCTGTAGGTCAGGTTGCTTTCTCAGGTATAGTATGAGGTCGCTCTTATTGACCACACC(.*)","\\1",df_PRA2,perl=T)



pos_mut2<-pos_mut %>%
group_by(pos,wt,mut) %>%
dplyr::mutate(Freq=sum(pos)/pos)%>%unique()


ggplot(pos_mut2,aes(x=pos,y=wt,fill=mut)) +
  geom_tile()+xlim(c(0,10))









ggplot(a,aes(x=pos,y=num,color = mut))+
  geom_line()+
  facet_wrap(~mut)


aa <- a %>%
  mutate(suibian = ifelse(num == 0, 0, 1))

ggplot(aa,aes(x=pos,y=mut,fill=suibian)) +
geom_tile() +
scale_fill_gradient(low="red",high="yellow")


for (i in 1:nrow(pos_mut2)){
  if(pos_mut2[i,1] < 150){
  pos_mut2[i,5] <-"0-150"
}else if(pos_mut2[i,1]<300){
  pos_mut2[i,5] <-"150-300"
}else if(pos_mut2[i,1]<450){
  pos_mut2[i,5] <-"300-450"
}else if(pos_mut2[i,1]<600){
  pos_mut2[i,5] <-"450-600"
}else if(pos_mut2[i,1]<750){
  pos_mut2[i,5] <-"600-750"
}else{
  pos_mut2[i,5] <-"750-900"
}}

ggplot(pos_mut2,aes(x=pos,y=Freq,color = mut))+ geom_point(size=0.5)+facet_grid( mut ~ part,scale="free_x")





cc_sex <- pra2 %>%
  melt(id.vars=c("gene","chr")) %>%
  mutate(sex = ifelse(value == 0 ,  "female","male"))




style.print <- function() {
  theme_classic() +
    theme(legend.text=element_text(size=unit(12,"bigpts")),#图例的内容
          legend.key.size=unit(12,"bigpts"),
          legend.key=element_blank(),
          legend.title=element_text(size=unit(12,"bigpts"),face="bold"),
          axis.title=element_text(size=unit(14,"bigpts")),#横纵坐标的标题
          axis.text=element_text(size=unit(12,"bigpts")))#横纵坐标刻度
}


#stickness

pra <- aa_hyb_RSA0.2_surfaced %>%
   dplyr::group_by(aa_hyb_RSA0.2_surfaced$AA_pos) %>%
   dplyr::summarise(hydrophobicity.y = mean(hydrophobicity.y))


hyb_RSA0.2_6_17 <- aa_hyb_RSA0.2_surfaced

# hyb_RSA0.2_6_17<- mutate(hyb_RSA0.2_6_17, 
#            if(hydrophobicity.y >0 & hydrophobicity.x <0){stickness = 0.5706377} 
#            else if(hydrophobicity.y <0 & hydrophobicity.x >0){stickness = 0.5501759}
#            else{stickness = 0.5683153})


hyb_RSA0.2_6_17 <-hyb_RSA0.2_6_17%>% mutate(stickness = ifelse(hydrophobicity.y >0 & hydrophobicity.x <0, 0.1142755, 
                                             ifelse(hydrophobicity.y <0 & hydrophobicity.x >0,0.108238, 0.111091)))

hyb_RSA0.2_6_17 <- hyb_RSA0.2_6_17 %>%
  mutate(sti_RSA = stickness*RSA)

hyb_RSA0.2_6_17 <- hyb_RSA0.2_6_17 %>%
  mutate(sti_ACC = stickness*ACC)

ggplot(hyb_RSA0.2_6_17,aes(hyb_RSA0.2_6_17$sti_RSA,hyb_RSA0.2_6_17$fitness)) +
  geom_point() +
  geom_smooth() + xlab("ch") + ylab("fitness")+ style.print()

cor.test(hyb_RSA0.2_6_17$sti_RSA,hyb_RSA0.2_6_17$fitness,method = "s")

ggplot(hyb_RSA0.2_6_17,aes(hyb_RSA0.2_6_17$sti_ACC,hyb_RSA0.2_6_17$fitness)) +
  geom_point() +
  geom_smooth() + xlab("chan") + ylab("fitness")+ style.print()

cor.test(hyb_RSA0.2_6_17$sti_ACC,hyb_RSA0.2_6_17$fitness,method = "s")



aa14 <- 0:10
aaa14 <- as.data.frame(aa14)

aaa14 <- aaa14 %>%
  mutate(p = dbinom(aa14,10,0.14))

ggplot(aaa14,aes(x=aa14,y=p)) +
   geom_point() +
  geom_smooth()

model <- loess(aaa14$p~aaa14$aa14)
x0 <- 1.4
y0 <- predict(model,x0)

#15
aa15 <- 0:10
aaa15 <- as.data.frame(aa15)

aaa15 <- aaa15 %>%
  mutate(p = dbinom(aa15,10,0.15))

ggplot(aaa15,aes(x=aa15,y=p)) +
  geom_point() +
  geom_smooth()

model <- loess(aaa15$p~aaa15$aa15)
x1 <- 1.5
y1 <- predict(model,x1)

#16
aa16 <- 0:10
aaa16 <- as.data.frame(aa16)

aaa16 <- aaa16 %>%
  mutate(p = dbinom(aa16,10,0.16))

ggplot(aaa16,aes(x=aa16,y=p)) +
  geom_point() +
  geom_smooth()

model <- loess(aaa16$p~aaa16$aa16)
x2 <- 1.6
y2 <- predict(model,x2)



a <- 0:5
aa <- as.data.frame(a)

aa <- aa %>%
  mutate(p = dbinom(a,5,0.14))

ggplot(aa,aes(x=a,y=p)) +
  geom_point() +
  geom_smooth()

# 处理数据

li <- seq(from=0.25, to=0.95, by=0.05)
lii <- as.data.frame(li)

ratio <- c(13/21,13/20,13/19,13/19,11/17,11/17,11/16,9/14,9/13,8/12,8/12,6/9,6/8,6/8,6/8)
ratio_df <- as.data.frame(ratio)

KC <- c(13,13,13,13,11,11,11,9,9,8,8,6,6,6,6)
KC_df <- as.data.frame(KC)

data_filt_df <- cbind(lii,ratio_df,KC_df)
names(data_filt_df) <- c("cutoff","ratio","KC")


ggplot(data_filt_df,aes(x=KC, y=ratio)) +
  geom_point() +
  geom_smooth(se= FALSE) + xlab("KC") + ylab("ratio")+ style.print()
  
ggplot(data_filt_df,aes(x=cutoff, y=ratio)) +
  geom_point() +
  geom_smooth(se= FALSE) + xlab("cutoff") + ylab("ratio")+ style.print()


