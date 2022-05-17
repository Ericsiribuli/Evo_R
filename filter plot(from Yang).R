library(plyr);
library(dplyr);
## number of mutant
dfBlast <- read.delim("parsed0_24.txt",header=F,sep="\t",stringsAsFactors = F);
dfBlast$varReg <- dfBlast$V3 %>%
  strsplit("")%>% 
  lapply(function(x){paste(x[15:83],collapse="",sep="")}) %>%
  unlist();
dfBlast$nMut <- dfBlast$V2 %>%
  strsplit("") %>% 
  lapply(function(x){length(which(x[15:83]==" "))}) %>%
  unlist();
dfBlast2 <- dfBlast %>%
  group_by(varReg) %>%
  sample_n(1);

dfSeqCount <- read.delim("dfT0_24.csv",sep=",",header=T,stringsAsFactors = F);
dfSeqCount$varReg <- 
  dfSeqCount$seq %>%
  strsplit("") %>% 
  lapply(function(x){paste(x[15:83],collapse="",sep="")}) %>%
  unlist();
dfSeqCount$fiveFixed <- 
  dfSeqCount$seq %>%
  strsplit("") %>% 
  lapply(function(x){paste(x[1:14],collapse="",sep="")}) %>%
  unlist();
dfSeqCount$threeFixed <- 
  dfSeqCount$seq %>%
  strsplit("") %>% 
  lapply(function(x){paste(x[84:89],collapse="",sep="")}) %>%
  unlist();
dfSeqCount2 <- dfSeqCount %>%
  filter(fiveFixed == "TTCAACCAAGTTGG") %>%
  filter(threeFixed == "CGTTGA") %>%
  group_by(varReg) %>%
  dplyr::summarise(Freq0 = sum(Freq0),
                   Freq24 = sum(Freq24));

dfSummary <- dfSeqCount2 %>% merge(dfBlast2,by.x='varReg',by.y='varReg')
dfSingle <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 1) %>%
  filter(!grepl("-",varReg));
dfWt <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 0) %>%
  filter(!grepl("-",varReg));

#double
dfDouble <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 2) %>%
  filter(!grepl("-",varReg));
#3-10mut
df3 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 3) %>%
  filter(!grepl("-",varReg))
df4 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 4) %>%
  filter(!grepl("-",varReg))
df5 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 5) %>%
  filter(!grepl("-",varReg))
df6 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 6) %>%
  filter(!grepl("-",varReg))
df7 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 7) %>%
  filter(!grepl("-",varReg))
df8 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut == 8) %>%
  filter(!grepl("-",varReg))
df9 <- dfSummary %>%
  filter(Freq0 >= 100) %>%
  filter(nMut >=9) %>%
  filter(!grepl("-",varReg))


dfSingle <- dfSingle %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
#double
dfDouble <- dfDouble %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-dfDouble$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
colnames(te)<-c("fitness")
dfDouble<-dfDouble[,-8]
dfDouble<-data.frame(dfDouble,te)


#df2 fig2E
te<-dfDouble$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df2<-data.frame(dfDouble,te)

#3-10mut
df3 <- df3 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df3$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df3<-data.frame(df3,te)

df4 <- df4 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df4$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df4<-data.frame(df4,te)

df5 <- df5 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df5$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df5<-data.frame(df5,te)

df6 <- df6 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df6$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df6<-data.frame(df6,te)


df7 <- df7 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df7$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df7<-data.frame(df7,te)

df8 <- df8 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df8$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df8<-data.frame(df8,te)

df9 <- df9 %>%
  mutate(fitness = (Freq24 / Freq0 / dfWt$Freq24 * dfWt$Freq0) ** (1/11.5) );
te<-df9$fitness
te<-ifelse(te<0.5,print(0.5),print(te))
te<-data.frame(te)
df9<-data.frame(df9,te)




dfSingle$mutPos <- dfSingle$V2 %>%
  strsplit("") %>% 
  lapply(function(x){which(x[15:83]==" ")}) %>%
  unlist();

dfDouble$mutPos <- dfDouble$V2 %>%
  strsplit("") %>% 
  lapply(function(x){which(x[15:83]==" ")}) %>% as.character()%>%
  unlist();

#double
SS <- dfDouble$V2 %>%
      strsplit("") %>% 
lapply(function(x){which(x[15:83]==" ")}) %>%
      unlist();



dfSingle <- dfSingle %>%
  group_by(varReg) %>%
  mutate(mutInto = strsplit(varReg,"")[[1]][mutPos])



library(ggplot2);
dfSingle %>%
  mutate(mutInto = factor(mutInto,levels=c("A","T","C","G"),ordered=T)) %>%
  ggplot(aes(x=mutPos,y=mutInto,fill=fitness)) +
  geom_tile() +
  scale_fill_gradient(low="red",high="yellow")

#double
ggplot(data=dfDoubleLa4,aes(x=a,y=b,fill=ε)) +
  geom_tile() +
  scale_fill_gradient(low="red",high="yellow")

