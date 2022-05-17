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

# wt_total <- data.frame(freq0=c("1654198"),freq7=c("2104566"))

wt_summary <- wt_summary %>%
  mutate(fitness = (freq7 / freq0 / 1230107 * 936831) ** (1/11.5) )

wt_summary <- wt_summary %>%
  filter(freq0>0)



wt_summary50 <- wt_summary %>%
  filter(freq0>50)
wt_summary100 <- wt_summary %>%
  filter(freq0>100)
wt_summary50$id <- row.names(wt_summary50)
wt_summary100$id <- row.names(wt_summary100)

ggplot(data=wt_summary100, aes(x=id,y=fitness)) +geom_point()
ggplot(data=wt_summary50, aes(x=id,y=fitness)) +geom_point()

wt_freq0_100_meansd <- wt_summary100%>%
  filter(fitness>0.83,fitness<1.02)
wt_freq0_100_meansd$id <- row.names(wt_freq0_100_meansd)
ggplot(data=wt_freq0_100_meansd, aes(x=id,y=fitness)) +geom_point() +ylim(0,1.5)

wt_freq0_50_meansd <- wt_summary50%>%
  filter(fitness>0.79,fitness<1.04)
wt_freq0_50_meansd$id <- row.names(wt_freq0_50_meansd)
ggplot(data=wt_freq0_50_meansd, aes(x=id,y=fitness)) +geom_point() +ylim(0,1.5)

wt_summary200 <- wt_summary %>%
  filter(freq0>200)
wt_summary200$id <- row.names(wt_summary200)


#single

s1_single <- read.csv("S1_single.csv",header = F)
names(s1_single) <- c("umi","pos")
s1_single$freq0 <- 1

s1_single <- s1_single %>%
  group_by(pos) %>%
  dplyr::summarise(freq0=sum(freq0))

s3_single <- read.csv("S3_single.csv",header = F)
names(s3_single) <- c("umi","pos")
s3_single$freq7 <- 1

s3_single <- s3_single %>%
  group_by(pos) %>%
  dplyr::summarise(freq7=sum(freq7))

single_summary <- merge(s1_single,s3_single,by.x="pos",all = T)
single_summary[is.na(single_summary)] <-0
single_summary1 <- single_summary %>%
  filter(freq0>100)

single_summary1 <- single_summary1 %>%
  mutate(fitness = (freq7 / freq0 / 661055 * 1071864) ** (1/11.5) )
single_summary1$id <- row.names(single_summary1)
ggplot(data=single_summary1, aes(x=id,y=fitness)) +geom_point() +ylim(0,1.5)

#3-25
s1_single2 <- s1_single %>%
  dplyr::group_by(umi) %>%
  dplyr::summarise(freq0=sum(freq0),pos=unique(pos))

s3_single2 <- s3_single %>%
  dplyr::group_by(umi) %>%
  dplyr::summarise(freq7=sum(freq7),pos=unique(pos))

single_summary <- merge(s1_single2,s3_single2,by.x="umi",by.y="umi",all = T)

single_summary[is.na(single_summary)] <-0
single_summary <- single_summary %>%
  filter(freq0>0)

single_summary1<-single_summary%>%
  mutate(fitness = (freq7 / freq0 / 1230107 * 936831) ** (1/11.5))


#correlation3-26
tx<-single_summary1

lapply(1:10, function(i){
  myData<- filter(tx,pos.x==tx$pos.x[i])
  myData<-myData[order(myData$freq0),]
  for (k in 1:nrow(myData)) {
    
    test0<-sum(myData$freq0[1:k])
    if(test0>=sum(myData$freq0)/2){
      next;
    }
    
    new1<-data.frame(sum0=test0,Time=k)
    
  }
  
  Time<-new1$Time+1
  TopF0<-filter(myData,myData$freq0[Time]>=freq0)
  BOTF0<-setdiff(myData,TopF0)
  
  CO_TOP<-cor.test(TopF0$freq0,BOTF0$freq0)
  CO_BOT<-cor.test(TopF0$freq7,BOTF0$freq7)
  fi<-data.frame(COTOP=CO_TOP$estimate,COBOT=CO_BOT$estimate)
})%>%rbind.fill()


#

s1_single3 <- s1_single %>%
  dplyr::group_by(pos) %>%
  dplyr::summarise(freq0 = sum(freq0))

s3_single3 <- s3_single %>%
  dplyr::group_by(pos) %>%
  dplyr::summarise(freq7 = sum(freq7))

single_summary2 <- merge(s1_single3,s3_single3,by = "pos",all = T)




