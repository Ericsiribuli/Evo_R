ccs_singlepos_mut<-read.csv("2-18singele_mut_ccs.csv",header = F)
names(ccs_singlepos_mut)<-c("pos","wt","mut")

ccs_singlepos_mut2<-ccs_singlepos_mut %>%
group_by(pos,wt,mut) %>%
dplyr::mutate(Freq_ccs=sum(pos)/pos)%>%unique()

ccs_singlepos_mut3 <- ccs_singlepos_mut2 %>%
  group_by(pos) %>%
  dplyr::summarise(Freq_ccs=sum(Freq_ccs))  


water_singlepos_mut<-read.csv("single_pos_mut.csv",header = F)
names(water_singlepos_mut)<-c("pos","wt","mut")

water_singlepos_mut2<-water_singlepos_mut %>%
  group_by(pos,wt,mut) %>%
  dplyr::mutate(Freq_water=sum(pos)/pos)%>%unique()

water_singlepos_mut3 <- water_singlepos_mut2 %>%
  group_by(pos) %>%
  dplyr::summarise(Freq_water=sum((Freq_water)) )

                    
singlepos_mut <- merge(ccs_singlepos_mut3,water_singlepos_mut3,all=T)
singlepos_mut[is.na(singlepos_mut)] <-0
singlepos_mut <- singlepos_mut%>%
  mutate(ccs_greater = Freq_ccs - Freq_water)

a <- which(singlepos_mut$ccs_greater>0)
length(a)
b <- which(singlepos_mut$ccs_greater<0)
length(b)
c <- which(singlepos_mut$ccs_greater==0)
length(c)

for (i in 1:nrow(singlepos_mut2)){
if(singlepos_mut2[i,1] < 160){
  singlepos_mut2[i,5] <-"0-160"
}else if(singlepos_mut2[i,1]<350){
  singlepos_mut2[i,5] <-"160-320"
}else if(singlepos_mut2[i,1]<480){
  singlepos_mut2[i,5] <-"320-480"
}else if(singlepos_mut2[i,1]<640){
  singlepos_mut2[i,5] <-"480-640"
}else if(singlepos_mut2[i,1]<800){
  singlepos_mut2[i,5] <-"640-800"
}else{
  singlepos_mut2[i,5] <-"800-960"
}}

names(ccs_singlepos_mut2)<-c("pos","wt","mut","Freq")
ggplot(ccs_singlepos_mut2,aes(x=pos,y=Freq,color = mut))+ geom_point(size=1)+facet_grid( mut ~.)+scale_x_continuous(breaks=seq(0, 1000, 50),expand = c(0, 0))

ggsave("2-18ccs_singlepos_point.png")


ggplot(ccs_singlepos_mut,aes(x=pos,color = mut))+ geom_bar(aes(fill = mut)) + scale_x_continuous(breaks=seq(0, 1000, 50),expand = c(0, 0))+ scale_y_continuous(breaks=seq(0, 90, 10),expand = c(0, 0))