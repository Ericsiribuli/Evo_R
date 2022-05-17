singlecell_data_Y <- filter(singlecell_data,chr=="Y")
singlecell_data_melt <- melt(singlecell_data,id.vars=c("gene","chr"))

tet <- singlecell_data_melt
tet$sex = rep(0,nrow(tet))
for (i in 1:nrow(tet)) {
  if(tet[i,4]==0){
    tet[i,5]<-"female"
  }else{
    tet[i,5]<-"male"
  }
}


CellName<-names(singlecell_data)[1:202]
FinalData<-lapply(3:ncol(singlecell_data),function(k){
  EachCell<-select(singlecell_data,c(1:2,k))
  names(EachCell)[3]<-"Exp"
  #计算有表达的X连锁基因的表达量中值Mx 
   MxData<-EachCell%>%
    filter(.,chr=="X")%>%
    mutate(ALLXGene=length(.$gene))%>%
    filter(.,Exp>0)%>%
    mutate(ExpXGene=length(.$gene))%>%
    summarise(cell=CellName[k],Mx=median(Exp),Fre=ExpXGene[1]/ALLXGene[1])
   #计算N%最高表达量的常染色体基因的表达量中值Ma 
   MaData<-EachCell%>%
     dplyr::filter(.,chr!="X" & chr!="Y")%>%
     arrange(.,desc(Exp))%>%
     dplyr::filter(.,Exp>=(quantile(.$Exp,1-MxData[,3])))%>%
     summarise(cell=CellName[k],Ma=median(Exp))
   FiData<-cbind(MxData,MaData)
   FiData$ratio=FiData$Mx/FiData$Ma
   return(FiData)
})%>%rbind.fill()



FinalData1 <- merge(FinalData,tet[,c(3,5)],by.x = "cell",by.y = "variable")

FinalData2 <- FinalData1
FinalData2$type <- rep(0,nrow(FinalData2))
FinalData2$time <- rep(0,nrow(FinalData2))

for (i in 1:nrow(FinalData2)) {
  FinalData2[i,7] <- unlist(strsplit(FinalData2[i,1],split = "[.]"))[2]
  FinalData2[i,8] <- unlist(strsplit(FinalData2[i,1],split = "[.]"))[[1]]
}

FinalData2$type <- as.factor(FinalData2$type)
FinalData2$time <- as.factor(FinalData2$time)

p1<-ggplot(FinalData2,aes(sex,ratio,color=sex)) +
    geom_boxplot()+theme_bw()
p2<-ggplot(FinalData2,aes(time,ratio,color=time)) +
    geom_boxplot()+theme_bw()
p3<-ggplot(FinalData2,aes(type,ratio,color=type)) +
    geom_boxplot()+theme_bw()
picture<-grid.arrange(p1,p2,p3,ncol=3)
 