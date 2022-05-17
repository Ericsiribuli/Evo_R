# 挑选TGY  day7&day5 fitness高和低各十个做实验-----

TGY_day5_day7 <- merge(TGY_day5,TGY_day7,by='mut',all=T)
names(TGY_day5_day7) <- c('mut','fitness_day5',"fitness_day7")

ggplot(TGY_day5_day7,aes(x =TGY_day5_day7$fitness_day5,TGY_day5_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGY_fitness_day5") + ylab("TGY_fitness_day7") +
  style.print()

TGY_day5_day3 <- merge(TGY_day5,TGY_day3,by='mut',all=T)
names(TGY_day5_day3) <- c('mut','fitness_day5',"fitness_day3")

ggplot(TGY_day5_day3,aes(x =TGY_day5_day3$fitness_day3,TGY_day5_day3$fitness_day5)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGY_fitness_day3") + ylab("TGY_fitness_day5") +
  style.print()

TGY_day3_day7 <- merge(TGY_day3,TGY_day7,by='mut',all=T)
names(TGY_day3_day7) <- c('mut','fitness_day3',"fitness_day7")

ggplot(TGY_day3_day7,aes(x =TGY_day3_day7$fitness_day3,TGY_day3_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGY_fitness_day3") + ylab("TGY_fitness_day7") +
  style.print()


TGS_day5_day7 <- merge(TGS_day5,TGS_day7,by='mut',all=T)
names(TGS_day5_day7) <- c('mut','fitness_day5',"fitness_day7")
TGS_day5_day7_te <- filter(TGS_day5_day7,fitness_day5<1&fitness_day7<1)

ggplot(TGS_day5_day7,aes(x =TGS_day5_day7$fitness_day5,TGS_day5_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.95,1.1) + ylim(0.95,1.1) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("TGS_fitness_day5") + ylab("TGS_fitness_day7") +
  style.print()

TGY_day5_day7_te <- filter(TGY_day5_day7,fitness_day5<0.996&fitness_day7<1)

TGY_TGS_day5_day7_te <- merge(TGY_day5_day7_te,TGS_day5_day7_te,by='mut',all = F)



AGY_day5_day7 <- merge(AGY_day5,AGY_day7,by='mut',all=T)
names(AGY_day5_day7) <- c('mut','fitness_day5',"fitness_day7")

ggplot(AGY_day5_day7,aes(x =AGY_day5_day7$fitness_day5,AGY_day5_day7$fitness_day7)) +
  geom_point(size=1) + xlim(0.98,1.03) + ylim(0.98,1.03) +
  geom_hline(aes(yintercept=1), linetype = "dashed") +
  geom_vline(aes(xintercept=1), linetype = "dashed") + geom_abline(intercept = 0, slope = 1 ) +xlab("AGY_fitness_day5") + ylab("AGY_fitness_day7") +
  style.print()



