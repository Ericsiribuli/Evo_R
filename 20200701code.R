Test<-readRDS("~/200331_fitness/AGY_geno.rds")
Test1<-Test%>%
  dplyr::filter(.,RSA >0.2 & hyb_t != 0 )


save.image("~/200331_fitness/7-1.Rdata")
load("~/200331_fitness/7-1.Rdata")