#200727 导入数据画出样本间相关性热图-----

library(ggplot2)
library(scales)


#TGY----
TGY_pearson <- read.csv("TGY_pearson.csv",header = F)
names(TGY_pearson) <- c("x","y","pearson")

TGY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PTDH3-GFP-SC Repeatability of freq",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#TGS----
TGS_pearson <- read.csv("TGS_pearson.csv",header = F)
names(TGS_pearson) <- c("x","y","pearson")

TGS_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PTDH3-GFP-SC Repeatability of freq",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

scale_x_discrete(breaks = seq(52,63,1),labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
scale_y_discrete(breaks = seq(52,63,1),labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))

  
#AGY----

AGY_pearson <- read.csv("AGY_pearson.csv",header = F)
names(AGY_pearson) <- c("x","y","pearson")

AGY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PADP1-GFP-YPD Repeatability of freq",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

  
#AGS ----
AGS_pearson <- read.csv("AGS_pearson.csv",header = F)
names(AGS_pearson) <- c("x","y","pearson")

AGS_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AGS_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#TUS-----
TUS_pearson <- read.csv("TUS_pearson.csv",header = F)
names(TUS_pearson) <- c("x","y","pearson")

TUS_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "TUS_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#TUU-----
TUU_pearson <- read.csv("TUU_pearson.csv",header = F)
names(TUU_pearson) <- c("x","y","pearson")

TUU_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "TUU_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#AUS-----
AUS_pearson <- read.csv("AUS_pearson.csv",header = F)
names(AUS_pearson) <- c("x","y","pearson")

AUS_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AUS_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#AUU-----
AUU_pearson <- read.csv("AUU_pearson.csv",header = F)
names(AUU_pearson) <- c("x","y","pearson")

AUU_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AUU_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#AUY-----
AUY_pearson <- read.csv("AUY_pearson.csv",header = F)
names(AUY_pearson) <- c("x","y","pearson")

AUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AUY_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#TUY-----
TUY_pearson <- read.csv("TUY_pearson.csv",header = F)
names(TUY_pearson) <- c("x","y","pearson")

p_TUY_pearson <-TUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "TUY_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


# nowt-------------

#TGY-----
TGY_pearson_nowt <- read.table("TGY_pearson_nowt.txt",sep = ',', stringsAsFactors = F)
names(TGY_pearson_nowt) <- c("x","y","pearson")

TGY_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PTDH3-GFP-YPD Repeatability of freq(No WT)",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#TGS-----

TGS_pearson_nowt <- read.csv("TGS_pearson_nowt.txt",sep = ',', stringsAsFactors = F,header = F)
names(TGS_pearson_nowt) <- c("x","y","pearson")

TGS_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "PTDH3-GFP-SC Repeatability of freq",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

scale_x_discrete(breaks = seq(52,63,1),labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(breaks = seq(52,63,1),labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))
  
#AGY -----

AGY_pearson_nowt <- read.csv("AGY_pearson_nowt.txt",sep = ',', stringsAsFactors = F,header = F)
names(AGY_pearson_nowt) <- c("x","y","pearson")

AGY_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AGY_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#AGS -----

AGS_pearson_nowt <- read.csv("AGS_pearson_nowt.txt",sep = ',', stringsAsFactors = F,header = F)
names(AGS_pearson_nowt) <- c("x","y","pearson")

AGS_pearson_nowt %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AGS_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#TUS-----
TUS_pearson <- read.csv("TUS_pearson.csv",header = F)
names(TUS_pearson) <- c("x","y","pearson")

TUS_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "TUS_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#TUU-----
TUU_pearson <- read.csv("TUU_pearson.csv",header = F)
names(TUU_pearson) <- c("x","y","pearson")

TUU_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "TUU_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#AUS-----
AUS_pearson <- read.csv("AUS_pearson.csv",header = F)
names(AUS_pearson) <- c("x","y","pearson")

AUS_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AUS_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')

#AUU-----
AUU_pearson <- read.csv("AUU_pearson.csv",header = F)
names(AUU_pearson) <- c("x","y","pearson")

AUU_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient(low="white",high="red") +
  labs(title = "AUU_pearson samples correlation",subtitle = "D1-D7 Timepoint,   S1-S3 Biological Repeat")+   ylab('') + xlab('')


#AUY----- 10-15version
AUY_pearson <- read.csv("AUY_pearson.csv",header = F)
names(AUY_pearson) <- c("x","y","pearson")

p_TUY_pearson <- TUY_pearson %>%
  mutate(y = factor(y,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  mutate(x = factor(x,levels=c("s1-d1","s2-d1","s3-d1","s1-d3","s2-d3","s3-d3","s1-d5","s2-d5","s3-d5","s1-d7","s2-d7","s3-d7"),ordered=T)) %>%
  ggplot(aes(x=x,y=y,fill=pearson)) +
  geom_tile() +
  scale_x_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_y_discrete(labels = c("D1-S1","D1-S2","D1-S3","D3-S1","D3-S2","D3-S3","D5-S1","D5-S2","D5-S3","D7-S1","D7-S2","D7-S3"))+
  scale_fill_gradient2(low="#FFE4E1",high="red",mid = '#FF7C74',midpoint = 0.994,limits=c(0.99,1),oob=squish,breaks=c(0.99,0.995,1), 
                      labels=c('<0.99',0.995,1)) + ylab('') + xlab('') +
  theme(axis.text.x= element_text(size=25,angle = 45,vjust = 0.5),axis.text.y= element_text(size=25),text= element_text(size=20))+theme(panel.border =element_blank())

p3 <- p_TUY_pearson +theme(panel.border =element_blank())

pdf("~/fig1I_remove_backg.pdf",width = 15,height = 14)
dev.off()



#TUY-----
TUY_pearson <- read.csv("TUY_pearson.csv",header = F)
names(TUY_pearson) <- c("x","y","pearson")

TUY_pearson %>%
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
                       labels=c(0.8,0.99,1),oob=squish)+ ylab('') + xlab('') +
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size = 18))+theme(panel.border =element_blank())

#######-------
scale_fill_gradientn(colours=c("white","#FFE4E1","#FFE4E1","red"),
                     values=c(0.8,0.9899,0.99,1),
                     na.value="white", guide="colourbar",
                     name="Pearson",limits=c(0.8,1),breaks=c(0.8,0.99,1), 
                     labels=c(0.8,0.99,1),oob=squish)


