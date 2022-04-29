# Load some packages 
install.packages("gamlss")
install.packages("emmeans")
install.packages("multcompView")
install.packages("moments")
install.packages("nortest")
install.packages("e1071")
install.packages("ggplot2")
install.packages("tibble")
library(gamlss)
library(emmeans)
library(multcompView)
library(moments)
library(nortest)
library(e1071)
library(ggplot2)
library(tibble)

# Attach file and organize the data table
library(dplyr)
leaf_mort <- read.csv('Leaf_mortality.csv', sep = ',')
leaf_mort$pop <- as.factor(leaf_mort$pop)

ntrt <- filter(leaf_mort, pop %in% c("SUS_con", "Flub_con", "H1_con", "H2_con"))
ntrt$pop <- factor(ntrt$pop, levels = c("SUS_con", "Flub_con", "H1_con", "H2_con"))

trt <- filter(leaf_mort, pop %in% c("SUS_flub", "Flub_flub", "H1_flub", "H2_flub"))
trt$pop <- factor(trt$pop, levels = c("SUS_flub", "Flub_flub", "H1_flub", "H2_flub"))

################################################################################
### GAMLSS
mod.ntrt <- gamlss(cbind(resp, total-resp)~pop, family = BI, what = "mu",
                   data = ntrt)
summary(mod.ntrt)

hist(mod.ntrt$residuals, freq = F)
lines(x = density(x = mod.ntrt$residuals), col = "red")

#########################################################
### Identify the best model that fits the data.

obj <- chooseDist(mod.ntrt, k = c(2, 4, 6), type = "binom",
                  extra = NULL, trace = FALSE,
                  parallel = "multicore")
getOrder(obj)

# Make update of the model.

fm<-update(mod.ntrt, family=names(getOrder(obj)[1]))

summary(fm)

hist(fm$residuals, freq = F)
lines(x = density(x = fm$residuals), col = "red")

##################################################
# Normality test of  residuals from gamlss model, to insecticide treatments.
resi.ntrt <- mod.ntrt$residuals
qqnorm(resi.ntrt)

# skewness and kurtosis, they should be around (0,3)
skewness(resi.ntrt)
kurtosis(resi.ntrt)

# Shapiro-Wilks test
shapiro.test(resi.ntrt)

# Kolmogorov-Smirnov test
ks.test(resi.ntrt,"pnorm",mean(resi.ntrt),sqrt(var(resi.ntrt)))

# Anderson-Darling test
ad.test(resi.ntrt)

################################################################################
### Identify the mean, confidence intervals and perform treatment contrasts.
nor.ntrt <- emmeans(mod.ntrt, "pop", df = mod.ntrt$df.residual) 
summary(nor.ntrt)

# Check for differences between treatments
ntrt.join <- joint_tests(nor.ntrt) 
ntrt.join

# Show confidence interval of treatments and z-test
confit.ntrt <- confint(nor.ntrt, level = 0.95, type = "response")
print(confit.ntrt)

ntrt.contr <- pairs(regrid(nor.ntrt), reverse = T) # regrid(normal)
ntrt.contr.summ <- as.data.frame(cbind(summary(ntrt.contr)[,1], round(summary(ntrt.contr)[,2:6], 3)))
colnames(ntrt.contr.summ) <- c('Contrast', 'estimate', 'SE', 'df', 't.ratio', 'p-value')

write.csv(ntrt.contr.summ, file = "Summary_Leaf_mortality_Control.csv")

################################################################################
### GAMLSS
mod.trt <- gamlss(cbind(resp, total-resp)~pop, family = BI, what = "mu",
                   data = trt)
summary(mod.trt)

hist(mod.trt$residuals, freq = F)
lines(x = density(x = mod.trt$residuals), col = "red")

#########################################################
### Identify the best model that fits the data.

obj <- chooseDist(mod.trt, k = c(2, 4, 6), type = "binom",
                  extra = NULL, trace = FALSE,
                  parallel = "multicore")
getOrder(obj)

# Make update of the model.

fm<-update(mod.trt, family=names(getOrder(obj)[1]))

summary(fm)

hist(fm$residuals, freq = F)
lines(x = density(x = fm$residuals), col = "red")

##################################################
# Normality test of  residuals from gamlss model, to insecticide treatments.
resi.trt <- mod.trt$residuals
qqnorm(resi.trt)

# skewness and kurtosis, they should be around (0,3)
skewness(resi.trt)
kurtosis(resi.trt)

# Shapiro-Wilks test
shapiro.test(resi.trt)

# Kolmogorov-Smirnov test
ks.test(resi.trt,"pnorm",mean(resi.trt),sqrt(var(resi.trt)))

# Anderson-Darling test
ad.test(resi.trt)

################################################################################
### Identify the mean, confidence intervals and perform treatment contrasts.
nor.trt <- emmeans(mod.trt, "pop", df = mod.trt$df.residual) 
summary(nor.trt)

# Check for differences between treatments
trt.join <- joint_tests(nor.trt) 
trt.join

# Show confidence interval of treatments and z-test
confit.trt <- confint(nor.trt, level = 0.95, type = "unlink")
print(confit.trt)

trt.contr <- pairs(regrid(nor.trt), reverse = T) # regrid(normal)
trt.contr.summ <- as.data.frame(cbind(summary(trt.contr)[,1], round(summary(trt.contr)[,2:6], 3)))
colnames(trt.contr.summ) <- c('Contrast', 'estimate', 'SE', 'df', 't.ratio', 'p-value')

write.csv(trt.contr.summ, file = "Summary_Leaf_mortality_Insecticide.csv")

################################################################################
###### Barplot of Seasons
library(dplyr)
ntrt.g <- ntrt %>% group_by(pop) %>% summarise(resp = sum(resp), total = sum(total))
ntrt.g$prop <- ntrt.g$resp/ntrt.g$total*100
sus.sd <- 0.032*sqrt(10)*100
flub.sd <- 0.000153*sqrt(10)*100
h1.sd <- 0.000153*sqrt(10)*100
h2.sd <- 0.0327*sqrt(10)*100
ntrt.g$sd <- c(sus.sd, flub.sd, h1.sd, h2.sd)

ntrt.gg <- ggplot(ntrt.g) +
  geom_bar( aes(x=pop, y=prop), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=pop, ymin=prop-sd, ymax=prop+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  theme_classic() +
  #ggtitle("Untreated Leaves") +
  xlab("Untreated Leaves") +
  ylab("% Survival") +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100)) +
  scale_x_discrete(labels=c("SUS_con" = "SUS", "Flub_con" = "Flub-R","H1_con" = "H1", "H2_con" = "H2")) +
  theme(plot.title = element_text(size = 30, hjust =0.5, vjust = -5)) +
  theme(axis.title.x = element_text(size = 40)) +
  theme(axis.title.y = element_text(size = 40)) +
  theme(axis.text.x = element_text(size = 40)) +
  theme(axis.text.y = element_text(size = 40))

# https://www.cedricscherer.com/2019/08/05/a-ggplot2-tutorial-for-beautiful-plotting-in-r/#text

trt.g <- trt %>% group_by(pop) %>% summarise(resp = sum(resp), total = sum(total))
trt.g$prop <- trt.g$resp/trt.g$total*100
sus.sd <- 0.0565*sqrt(10)*10
flub.sd <- 0.0001537*sqrt(10)*100
h1.sd <- 0.03277306*sqrt(10)*100
h2.sd <- 0.00015373*sqrt(10)*100
trt.g$sd <- c(sus.sd, flub.sd, h1.sd, h2.sd)
trt.g$cld <- c("a", "b", "b", "b")

trt.gg <- ggplot(trt.g) +
  geom_bar( aes(x=pop, y=prop), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=pop, ymin=prop-sd, ymax=prop+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  theme_classic() +
  #ggtitle("Untreated Leaves") +
  xlab("Insecticide Treated Leaves") +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100)) +
  scale_x_discrete(labels=c("SUS_flub" = "SUS", "Flub_flub" = "Flub-R","H1_flub" = "H1", "H2_flub" = "H2")) +
  theme(plot.title = element_text(size = 30, hjust =0.5, vjust = -5)) +
  theme(axis.title.x = element_text(size = 40)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 40)) +
  theme(axis.text.y = element_text(size = 40)) 

install.packages("ggpubr")
library(ggpubr)

ggarrange(ntrt.gg, trt.gg,
          ncol = 2, nrow = 1, font.label=list(color="black",size=40)) # , labels = c("A", "B")
save.image(file = "Leaf_Mortality.RData")
