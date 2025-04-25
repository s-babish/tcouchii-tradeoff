#libraries ----
library(tidyverse)
library(lmtest)

#Runs stats on the x0 results from fitting sigmoidal curves to rheobase data
rheodat <- read.csv("OutFiles/Rheobase/Sigmoidal/test/Sigmoidal-rpt.csv")
rheodat$pulse_length <- as.factor(rheodat$pulse_length)

#anova on x0 values between pulses ----
rheo_anova <- aov(x0 ~ pulse_length, data = rheodat)
summary(rheo_anova)

rheo_tukey <- TukeyHSD(rheo_anova)
rheo_tukey
#no differences in x0 between any of the groups (at least not when combined)

# MAMU vs x0 ----
summary(glm(x0 ~ MAMU*pulse_length, data = rheodat))
#not related to either variable
#just trying to cover my bases so I don't get roasted for not doing the "right" stats

#forget pulse length and just look at the plain data
lm1 <- lm(x0 ~ MAMU, data = rheodat)
summary(lm1)
plot(lm1)
#still no effect, as expected

#Correlation test 
pearson_corr <- cor.test(rheodat$MAMU, rheodat$x0, method = "pearson")
pearson_corr
#not at all correlated

# MAMU vs change between pulse lengths (delta x0)
#drop the first row, pick just the columns we care about, arrange by snake and
# then group by it to calculate differences between pulse lengths (which are
# already in order)
deltadat <- rheodat[-1,] %>% 
  select(X,x0,MAMU,Snake,pulse_length) %>% 
  arrange(X) %>% 
  group_by(X) %>%
  mutate(deltax0 = x0 - lag(x0)) %>% 
  na.omit() #the ones with NAs don't tell us anything

deltadat

summary(glm(deltax0 ~ MAMU*pulse_length, data = deltadat))

rheo_anova2 <- aov(deltax0 ~ pulse_length, data = deltadat)
summary(rheo_anova2)

rheo_tukey2 <- TukeyHSD(rheo_anova2)
rheo_tukey2
#big difference in 5000-10000 and 10000-50000

#do some plotting to look at this
deltadat %>% ggplot( aes(x=pulse_length, y=deltax0, fill=pulse_length)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Rheobase delta x0 vs pulse width") +
  xlab("Pulse width (us)")

deltadat %>% 
  ggplot(aes(MAMU, deltax0, color = as.factor(pulse_length)))+
  geom_point()
#definitely nothing to do with MAMU

#*need to ask if the delta x0 is a real metric or something i invented and go from there ----

#above but with IC50 ----
#load data and combine with current data
IC50 <- read.csv("data_raw/IC50/IC50.csv")
IC50dat <- merge(rheodat,IC50, by="Snake")
IC50dat <- IC50dat[!is.na(IC50dat$IC50),]
IC50dat[IC50dat$x0 == 600,]
#Snake 3064, CRF3064_m1_10, pulse length 50000
#remove that guy for now
IC50dat <- IC50dat[IC50dat$x0 < 600,]
IC50dat$Snake <- as.factor(IC50dat$Snake)

unique(IC50dat[!is.na(IC50dat$IC50), 1])
#only 14 individuals in this group

#do a glm 
IC50_glm <- glm(x0 ~ IC50*pulse_length, data = IC50dat)
plot(IC50_glm)
summary(IC50_glm)

#residuals are heteroskedastic (but not with the 600 outlier removed)
bptest(IC50_glm)

#forget pulse length and just look at the plain data
lm2 <- lm(x0 ~ IC50, data = IC50dat)
summary(lm2)
plot(lm2)
#vague relationship between x0 and IC50 (which keeps appearing and disappearing?)
#Correlation test
pearson_corr2 <- cor.test(IC50dat$IC50, IC50dat$x0, method = "pearson")
pearson_corr2
#allegedly they are correlated???
#also it generally looks like they're negatively correlated, so higher IC50 = 
# lower x0 = more excitability (which is opposite what bobby suggested could happen)
#kierstin says they could both be being influenced by a secret third thing

#delta x0 data
IC50dat <- deltadat %>% 
  mutate( 
    IC50 = 0) %>% 
  rows_update(distinct(IC50[,c(2,10)]), by = "Snake", unmatched = "ignore") %>% 
  mutate (
    IC50 = ifelse(IC50 == 0, NA, IC50)
  )

summary(glm(deltax0 ~ IC50*pulse_length, data = IC50dat))
#very small effect of IC50 and longest pulse length on deltax0

IC50dat %>% 
  ggplot(aes(IC50, deltax0, color = as.factor(pulse_length)))+
  geom_point()
#I suspect that one high value of delta x0 is what's causing the relationships
# really need to see if that's a real parameter of interest or not