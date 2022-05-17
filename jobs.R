#3-1
pra <- read.csv("singlecell_data.txt",sep = "\t")


#生信编程作业选做3/2

pra2 <- pra %>% filter(chr == "Y")
new_pra <- pra2 %>%
  melt(id.vars=c("gene","chr")) %>%
  mutate(sex = ifelse(value == 0,"female","male"))



#生信编程作业选做3/3/1
Mx_x <- pra %>%
  filter(chr == "X") %>%
  melt(id.vars=c("gene","chr")) %>%
  filter(value != 0) %>%
  mutate(gene = "expression") %>%
  dcast(variable ~ gene, median, value.var = "value")


#生信编程作业选做3/3/2
x_express_gene_number <- pra %>%
  filter(chr == "X") %>% nrow()
X_list <- pra %>%
  filter(chr == "X") 
data_x <- X_list[,3:ncol(X_list)]
data_x_matrix <- as.matrix(data_x)
result_values <- lapply(1:ncol(data_x_matrix),function(x){
  ((x_express_gene_number-table(data_x_matrix[,x])[1])/x_express_gene_number)*100
})%>%unlist()
answer2 <- data.frame(Cell_names=colnames(data_x),Percentage=result_values)


#生信编程作业选做3/3/3
options(digits =4)
autochr_nums <- autochr <- aa %>%
  filter(chr != "X" & chr !="Y" )%>%
  nrow()
autochr <- pra %>%
  filter(chr != "X" & chr !="Y" )
autochr <- autochr[,3:ncol(autochr)]
result_values_auto<- lapply(1:ncol(autochr),function(y){
  threshold <- (answer2[answer2$Cell_names==(colnames(autochr)[y]),2]/100)*autochr_nums%>%
    ceiling()
  sort(autochr[,y],decreasing = T)[1:threshold]%>%median()
})%>%unlist()
answer3.1 <- data.frame(Cell_names=colnames(autochr),Medain_auto=result_values_auto)
answer3.2 <- merge(x_Mx,answer3.1,by.x ="variable",by.y = "Cell_names")%>%
  mutate(Ratio_Mx_Ma=expression/Medain_auto)


#生信编程作业选做3/4
for_draw <- answer3.2%>%
  mutate(phase=gsub("(E\\d+).(\\w+).*","\\1",variable),cell_type=gsub("(E\\d+).(\\w+).*","\\2",variable)) %>%
  merge(.,cc_sex[,c(3,5)],by = "variable")
boxplot(for_draw[,4]~for_draw[,5],outlier=F)
boxplot(for_draw[,4]~for_draw[,6],outlier=F)
boxplot(for_draw[,4]~for_draw[,7],outlier=F)


