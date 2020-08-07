library(survival) # for survival plots
library(tidyverse) # for wrangling

setwd("~/Dropbox/tsa_2014/survival/")

# load data and format for analysis
d <- read.csv("fgb_survival.csv")
head(d)
d$bank <- as.factor(d$bank)
d$geno <- as.factor(d$geno)
d$rep <- as.factor(d$rep)
summary(d)
str(d)

# replace clones with new ids
d$oldgeno <- d$geno
d$geno <- gsub("^7","729",d$geno)
d$geno <- gsub("^29","729",d$geno)
d$geno <- gsub("^21","2117",d$geno)
d$geno <- gsub("^17","2117",d$geno)

d$geno <- as.factor(d$geno)

d %>% filter(treat=="c") %>% nrow()
d %>% filter(treat=="c" & dead=="1") %>% nrow()
(3/56)*100
d %>% filter(treat=="v") %>% nrow()
d %>% filter(treat=="v" & dead=="1") %>% nrow()
(41/53)*100
# total dead = 
(3+41)/(56+53)
# 40 %

d %>% filter(bank=="e") %>% nrow()
d %>% filter(bank=="e" & dead=="1") %>% nrow()
(19/56)*100
d %>% filter(bank=="w") %>% nrow()
d %>% filter(bank=="w" & dead=="1") %>% nrow()
(25/53)*100
# total dead = 
(19+25)/(56+53)
# 40 %

# Survival by Treatment ------
die.treat <- Surv(d$tod, d$dead)
head(die.treat)
surT <- survfit(die.treat ~ treat, data=d)
summary(surT)

phTreat <- coxph(die.treat ~ bank + geno + treat, data=d)
summary(phTreat)
pVals <- anova(phTreat)

pBank_legend <- paste("p = ", round(pVals$`Pr(>|Chi|)`[2],2), sep="")
pGenet_legend <- paste("p = ", round(pVals$`Pr(>|Chi|)`[3],2), sep="")
pTreat_legend <- paste("p = ", round(pVals$`Pr(>|Chi|)`[4],2), sep="")

#       loglik   Chisq Df Pr(>|Chi|)    
# NULL  -196.28                          
# bank  -195.14  2.2848  1     0.1306    
# geno  -159.74 70.7867 25  2.938e-06 ***
# treat -120.14 79.2075  1  < 2.2e-16 ***

quartz()
par(mfrow=c(1,2))
plot(surT, col=c("grey60","red4"), lwd=5, xlab="Days", ylab="Fraction Surviving", main="Survival by Treatment")
abline(v=6,col="grey60", lty=2)
abline(v=10,col="grey60", lty=3)
legend("bottomleft", col=c("grey60","red4"), legend=c("Control", "Vibrio"), lwd=5, bty = 'n')
legend("bottomright", legend=pTreat_legend, bty = "n")

# Survival by bank ------
die.bank <- Surv(d$tod, d$dead)
head(die.bank)
surB <- survfit(die.bank ~ bank, data=d)
summary(surB)

cphBank <- coxph(die.bank ~ bank, data=d)

plot(surB, col=c("navy", "darkgreen"), lwd=5, xlab="Days", ylab="Fraction Surviving", main="Survival by Collection Site")
abline(v=6,col="grey60", lty=2)
abline(v=10,col="grey60", lty=3)
legend("bottomleft", col=c("navy", "darkgreen"), legend=c("East Bank", "West Bank"), lwd=5, bty = "n")
legend("bottomright", legend=pBank_legend, bty = "n")

# Make combined treatment + bank

head(d)
d <- d %>% mutate(bank_treat = paste(bank,treat,sep="_"))

die.bank_treat <- Surv(d$tod, d$dead)
head(die.bank_treat)
surBT <- survfit(die.treat ~ bank_treat, data=d)
summary(surBT)

phTreat <- coxph(die.treat ~ bank + geno + treat, data=d)
summary(phTreat)
pVals <- anova(phTreat)

pBank_legend <- paste("p(bank) = ", round(pVals$`Pr(>|Chi|)`[2],2), sep="")
pGenet_legend <- paste("p = ", round(pVals$`Pr(>|Chi|)`[3],2), sep="")
pTreat_legend <- paste("p(treatment) = ", round(pVals$`Pr(>|Chi|)`[4],2), sep="")

#       loglik   Chisq Df Pr(>|Chi|)    
# NULL  -196.28                          
# bank  -195.14  2.2848  1     0.1306    
# geno  -159.74 70.7867 25  2.938e-06 ***
# treat -120.14 79.2075  1  < 2.2e-16 ***

quartz()
plot(surBT, col=c("navy","navy", "darkgreen", "darkgreen"), lty=c(1,3,1,3),
     lwd=5, xlab="Days", ylab="Fraction Surviving", main="Survival by Treatment")
abline(v=6,col="grey60", lty=2)
abline(v=10,col="grey60", lty=3)
legend("bottomleft", col=c("black","black","navy","darkgreen"),lty=c(1,3,1,1), 
       legend=c("Control", "Vibrio", "East FGB", "West FGB"), lwd=5, bty = 'n')
legend("bottomright", legend=pTreat_legend, bty = "n")
legend(-1,0.5, legend=pBank_legend, bty = "n")
