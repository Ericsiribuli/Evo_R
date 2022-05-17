dfSingleONE<- dfSingle %>%
  group_by(mutPos) %>%
  dplyr::summarise(fitness=mean(fitness)) 

dfSingleONE2<- dfSingleONE

dfDoubleONE<- dfDouble %>%
  group_by(mutPos) %>%
  dplyr::summarise(fitness=mean(fitness)) 




#useless
S <- dfDouble$V2 %>%
  strsplit("") %>% 
  lapply(function(x){which(x[15:83]==" ")}) %>% as.numeric() %>%
  unlist();

da<-lapply(1:3,function(i){
   
   d<-cbind(SS[i],SS[i+1])
   d<-as.data.frame(d)
   i=i+2
})%>% rbind.fill()

#continue
a<-1
for (i in 0:11444) {
  a[i]=SS[i*2+1]
  unlist(a)
}

b<-1
for (i in 1:11445) {
  b[i]=SS[i*2]
  unlist(b)
}
qqa<-as.data.frame(a)
#a,b different length
row<-c(1)
qqa<-rbind(row,qqa)

qqb<-as.data.frame(b)

ab<-data.frame(qqa,qqb)

dfDouble2<-data.frame(dfDouble,ab)


dfDouble3<- dfDouble2 %>%
  group_by(mutPos) %>%
  dplyr::summarise(fitness=mean(fitness)
                   ,a=mean(a),b=mean(b))  


names(dfDouble3)<-c("mutPos","Fab","a","b")
names(dfSingleONE)<-c("a","Fa")
names(dfSingleONE2)<-c("b","Fb")

dfDoublete<-merge(x=dfDouble3,y=dfSingleONE,by="a",all=TRUE)
dfDoubleLa<-merge(x=dfDoublete,y=dfSingleONE2,by="b",all=TRUE)
dfDoubleLa<-dfDoubleLa[-1,]
dfDoubleLa2<-mutate(dfDoubleLa,FaFb=Fa*Fb)

x<-dfDoubleLa2$FaFb
xx<-ifelse(x<0.5,print(0.5),print(x))
xx<-data.frame(xx)
dfDoubleLa3<-data.frame(dfDoubleLa2,xx)
dfDoubleLa4<-mutate(dfDoubleLa3,ε=Fab-xx)
dfDoubleLa4<-dfDoubleLa4[-2334,]

ggplot(data=dfDoubleLa4,aes(x=a,y=b,fill=ε)) +
  geom_tile() +
  scale_fill_gradient2(low="red",high="blue",mid="yellow",breaks=c(-0.4,-0.2,0,0.2,0.4))

ggplot(data=dfDoubleLa4,aes(x=a,y=b,fill=ε)) +
  geom_tile() +
  scale_fill_gradientn(colours = c(colorRampPalette(c("blue","skyblue")),colorRampPalette(c("yellow","red"))),breaks=c(-0.4,-0.2,0,0.2,0.4))



#fig 2E
a<-data.frame(sd(dfSingle$fitness),min(dfSingle$fitness),mean(dfSingle$fitness))
b<-data.frame(max(df2$te),min(df2$te),mean(df2$te))
c<-data.frame(max(df3$te),min(df3$te),mean(df3$te))
d<-data.frame(max(df4$te),min(df4$te),mean(df4$te))
e<-data.frame(max(df5$te),min(df5$te),mean(df5$te))
f<-data.frame(max(df6$te),min(df6$te),mean(df6$te))
g<-data.frame(max(df7$te),min(df7$te),mean(df7$te))
h<-data.frame(max(df8$te),min(df8$te),mean(df8$te))
i<-data.frame(max(df9$te),min(df9$te),mean(df9$te))
names(a)<-c("MaxFitness","MinFitness","MeanFitness")
names(b)<-c("MaxFitness","MinFitness","MeanFitness")
names(c)<-c("MaxFitness","MinFitness","MeanFitness")
names(d)<-c("MaxFitness","MinFitness","MeanFitness")
names(e)<-c("MaxFitness","MinFitness","MeanFitness")
names(f)<-c("MaxFitness","MinFitness","MeanFitness")
names(g)<-c("MaxFitness","MinFitness","MeanFitness")
names(h)<-c("MaxFitness","MinFitness","MeanFitness")
names(i)<-c("MaxFitness","MinFitness","MeanFitness")
fitnessmm<-rbind.data.frame(a,b,c,d,e,f,g,h,i)
te<-as.data.frame(1:9)
names(te)<-c("NumMut")
fitnessmm<-data.frame(fitnessmm,te)

a<-data.frame(sd(dfSingle$fitness),mean(dfSingle$fitness))
b<-data.frame(sd(df2$te),mean(df2$te))
c<-data.frame(sd(df3$te),mean(df3$te))
d<-data.frame(sd(df4$te),mean(df4$te))
e<-data.frame(sd(df5$te),mean(df5$te))
f<-data.frame(sd(df6$te),mean(df6$te))
g<-data.frame(sd(df7$te),mean(df7$te))
h<-data.frame(sd(df8$te),mean(df8$te))
i<-data.frame(sd(df9$te),mean(df9$te))
names(a)<-c("sd","MeanFitness")
names(b)<-c("sd","MeanFitness")
names(c)<-c("sd","MeanFitness")
names(d)<-c("sd","MeanFitness")
names(e)<-c("sd","MeanFitness")
names(f)<-c("sd","MeanFitness")
names(g)<-c("sd","MeanFitness")
names(h)<-c("sd","MeanFitness")
names(i)<-c("sd","MeanFitness")
fitness_sd_fit<-rbind.data.frame(a,b,c,d,e,f,g,h,i)
tee<-as.data.frame(1:9)
names(tee)<-c("NumMut")
fitness_sd_fit<-data.frame(fitness_sd_fit,tee)
te<-fitness_sd_fit$MeanFitness-fitness_sd_fit$sd
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
fitness_sd_fit$Minerrorbar<-te
fitness_sd_fit<-fitness_sd_fit[,-4]

ggplot()+
  geom_point(data=fitness_sd_fit,mapping = aes(x=NumMut,y=MeanFitness))+
  geom_errorbar(data=fitness_sd_fit,mapping=aes(x=NumMut, ymin=Minerrorbar, ymax=MeanFitness+sd),width=.1)+
  scale_x_continuous(breaks=seq(0:9))+
  scale_y_continuous(breaks=seq(0.5:1.3,by=0.1))

ggplot()+
  geom_point(data=fitnessmm,mapping = aes(x=NumMut,y=MeanFitness))+
  geom_error(data=fitnessmm,mapping=aes(x=NumMut, ymin=MinFitness, ymax=MaxFitness))+
  scale_x_continuous(breaks=seq(0:9))+
  scale_y_continuous(breaks=seq(0.5:1.3,by=0.1))
