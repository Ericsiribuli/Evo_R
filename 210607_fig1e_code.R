P_Fig1E <- TGY_pearson_nowt_037 %>%
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
ss