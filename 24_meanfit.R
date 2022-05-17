library(ShortRead)
q948<-readFastq("SRR3159948.merged.fastq")
q948<-sread(q948)
q948<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q948,perl=T)
q948<-sub("(.*CGTTGATTA).*","\\1",q948,perl=T)
seq948<-as.data.frame(table(q948))

q945<-readFastq("SRR3159945.merged.fastq")
q945<-sread(q945)
q945<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q945,perl=T)
q945<-sub("(.*CGTTGATTA).*","\\1",q945,perl=T)
seq945<-as.data.frame(table(q945))

q944<-readFastq("SRR3159944.merged.fastq")
q944<-sread(q944)
q944<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q944,perl=T)
q944<-sub("(.*CGTTGATTA).*","\\1",q944,perl=T)
seq944<-as.data.frame(table(q944))

q941<-readFastq("SRR3159941.merged.fastq")
q941<-sread(q941)
q941<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q941,perl=T)
q941<-sub("(.*CGTTGATTA).*","\\1",q941,perl=T)
seq941<-as.data.frame(table(q941))

q868<-readFastq("SRR3159868.merged.fastq")
q868<-sread(q868)
q868<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q868,perl=T)
q868<-sub("(.*CGTTGATTA).*","\\1",q868,perl=T)
seq868<-as.data.frame(table(q868))

q840<-readFastq("SRR3159840.merged.fastq")
q840<-sread(q840)
q840<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q840,perl=T)
q840<-sub("(.*CGTTGATTA).*","\\1",q840,perl=T)
seq840<-as.data.frame(table(q840))

colnames(seq840)=c("seq","Freq840")
colnames(seq868)=c("seq","Freq868")
colnames(seq941)=c("seq","Freq941")
colnames(seq944)=c("seq","Freq944")
colnames(seq945)=c("seq","Freq945")
colnames(seq948)=c("seq","Freq948")
dfT24<-merge(seq840,seq868,all=T)
dfT24<-merge(dfT24,seq941,all=T)
dfT24<-merge(dfT24,seq944,all=T)
dfT24<-merge(dfT24,seq945,all=T)
dfT24<-merge(dfT24,seq948,all=T)
dfT24[is.na(dfT24)]<-0

#slowly
#me<-mclapply(1:nrow(dfT24),function(i){
#  tt<-mean(as.numeric(te[i,]))
# new<-data.frame(mean=tt,id=i)
#},mc.cores = 10)%>%rbind.fill()

te<-dfT24[,2:7]
tee<-data.frame(apply(te,1,mean))
dfT24$MeanFit24<-tee


#T0
q722<-readFastq("SRR3159722.merged.fastq")
q722<-sread(q722)
q722<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q722,perl=T)
q722<-sub("(.*CGTTGATTA).*","\\1",q722,perl=T)
seq722<-as.data.frame(table(q722))

q828<-readFastq("SRR3158828.merged.fastq")
q828<-sread(q828)
q828<-sub(".*(TTCAACCAAGTTGG.*)","\\1",q828,perl=T)
q828<-sub("(.*CGTTGATTA).*","\\1",q828,perl=T)
seq828<-as.data.frame(table(q828))

colnames(seq722)=c("seq","Freq722")
colnames(seq828)=c("seq","Freq828")
dfT0<-merge(seq722,seq828,all=T)
dfT0[is.na(dfT0)]<-0

te<-dfT0[,2:3]
tee<-data.frame(apply(te,1,mean))
dfT0$MeanFit0<-tee
dfT0<-subset(dfT0,MeanFit0>=100)
colnames(dfT0)=c("seq","Freq722","Freq828","FitnessT0")

dfT0_24<-merge(dfT24,dfT0,all.x = F,all.y = T)
dfte<-dfT0_24[,-2]
te[is.na(te)]<-0
#delete colo
colnames(te)=c("seq","Freq0","Freq24")
write.csv(te,"df0_24.csv")