#libraries ----
library(tidyverse)
library(lmtest)
library(minpack.lm)

foldername <- "WT_elegans"
#Runs stats on the x0 results from fitting sigmoidal curves to rheobase data
rheodat <- read.csv(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv",sep=""))
rheodat$pulse_length <- as.factor(rheodat$pulse_length)
rheodat <- rheodat[-1,-1]

#anova on x0 values between pulses ----
rheo_anova <- aov(x0 ~ pulse_length, data = rheodat)
summary(rheo_anova)

rheo_tukey <- TukeyHSD(rheo_anova)
rheo_tukey
#10000 is different from some of the others (100 and 50000)

#EC10 ----
ec10_anova <- aov(EC10 ~ pulse_length, data = rheodat)
summary(ec10_anova)

ec10_tukey <- TukeyHSD(ec10_anova)
ec10_tukey
#all the same 

#EC90 ----
#remove overly large ec90 values
rheodat <- rheodat %>% 
  filter(EC90<1000)

ec90_anova <- aov(EC90 ~ pulse_length, rheodat)
summary(ec90_anova)

ec90_tukey <- TukeyHSD(ec90_anova)
ec90_tukey
#lots of differences here but not with any clear trend so i won't trust it until
# cleaning issues in the raw data

#(but for now, it says 10000 is different from the others, the others are all identical )
#max force between pulses ----
rheo_anova2 <- aov(maxF ~ pulse_length, data = rheodat)
summary(rheo_anova2)

rheo_tukey2 <- TukeyHSD(rheo_anova2)
rheo_tukey2
#10000 and 50000 absolutely different from everything, and from each other

# MAMU vs x0 and max force ----
summary(glm(x0 ~ MAMU*pulse_length, data = rheodat))
#not related to either variable
#just trying to cover my bases so I don't get roasted for not doing the "right" stats

summary(glm(maxF ~ MAMU*pulse_length, data = rheodat))
#only related to pulse length, not MAMU or any interaction
summary(glm(maxF ~ MAMU, data = rheodat))
#slightly more suggestive when there's no pulse length splitting,
# but not nearly significant still

#forget pulse length and just look at the plain data
lm1 <- lm(x0 ~ MAMU, data = rheodat)
summary(lm1)
plot(lm1)
#still no effect, as expected

#Correlation test 
pearson_corr <- cor.test(rheodat$MAMU, rheodat$x0, method = "pearson")
pearson_corr
#not at all correlated

force_corr <- cor.test(rheodat$MAMU, rheodat$maxF, method = "pearson")
force_corr
#not correlated either

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

IC50_force <- glm(maxF ~ IC50*pulse_length, data = IC50dat)
summary(IC50_force)
plot(IC50_force)
#not great residuals for this but I'll worry about it later

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

#exponential decay curves fit to x0 (most current analysis) ----
all_fits <- matrix(ncol=4)
colnames(all_fits) <- c("Ind","x0","t","Genotype")

genotypes <- c("LVNV","WT_sirtalis","EPN","P","T","WT_elegans","WT_hammondii")

for (type in genotypes) {
  foldername = type
  df<-read.csv(toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv", sep="")))
  df <- df[-1,-1] %>% 
    filter(MusMassg != 0.000999) %>% 
    dplyr::select(X,x0,pulse_length) %>% 
    pivot_wider(names_from = pulse_length, values_from = x0) %>% 
    column_to_rownames('X') 
  df <- t(df)
  df <- df[1:5,] #ditching the high pulses bc they make fitting harder
  
  expdecay <- function(pw,x0,tau) {
    x0*exp(-pw*tau)
  }
  expdecay <- Vectorize(expdecay)
  
  expresidFun <- function(parS,observed,indices){
    expdecay(as.numeric(rownames(df))[indices],parS$x0,parS$tau) - observed
  }
  
  fitDecay <- function(x){
    nls.out <- nls.lm(par = list(x0 = df[1,x],tau = 0.1), fn = expresidFun, 
                      control = nls.lm.control(maxiter = 1024, maxfev = 1024), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
    print(nls.out)
    unlist(nls.out$par)
  }
  
  result <- Vectorize(fitDecay)(1:ncol(df))
  print(result)
  
  output <- t(matrix(c(colnames(df),result,rep(foldername,ncol(df))),nrow=4,byrow=T))
  colnames(output) <- c("Ind","x0","t","Genotype")
  all_fits <- rbind(all_fits,output)
}

all_fits <- as.data.frame(all_fits[-1,]) #not sure why it was formatted weirdly
#write.csv(all_fits,"OutFiles/Rheobase/exp_decay_method2.csv",row.names=F)
all_fits <- as.data.frame(all_fits)
comparison <- kruskal.test(all_fits$t ~ all_fits$Genotype)
print(comparison)

comparison2 <- kruskal.test(all_fits$x0 ~ all_fits$Genotype)
print(comparison2)

# dunn <- c("exp_decay_t",dunn.test::dunn.test(all_fits$t, g = all_fits$Genotype, 
#                                        kw = T, method = "bonferroni",
#                                        table = T, list = F))

dunn.test::dunn.test(as.numeric(all_fits$x0), g = all_fits$Genotype, 
          kw = T, method = "bonferroni",
          table = T, list = F)

#now do it again with the curve typically used for rheobase data -----
all_fits <- matrix(ncol=3)
colnames(all_fits) <- c("Ind","rheobase","Genotype")

genotypes <- c("LVNV","WT_sirtalis","T","WT_hammondii","EPN","WT_atratus")

for (type in genotypes) {
  foldername = type
  df<-read.csv(toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv", sep="")))
  df <- df[-1,-1] %>% 
    filter(MusMassg != 0.000999) %>% 
    dplyr::select(X,x0,pulse_length) %>% 
    pivot_wider(names_from = pulse_length, values_from = x0) %>% 
    column_to_rownames('X') 
  df <- t(df)
  df <- df[1:5,] #ditching the high pulses bc they make fitting harder
  
  expdecay <- function(pw,rheobase) {
    rheobase*(1+(2*rheobase/pw)) #chronaxie is just 2*rheobase
  }
  expdecay <- Vectorize(expdecay)
  
  expresidFun <- function(parS,observed,indices){
    expdecay(as.numeric(rownames(df))[indices],parS$rheobase) - observed
  }
  
  #parStart <- list(x0 = exp(coef(fm0)[[1]]), tau = -coef(fm0)[[2]])
  
  fitDecay <- function(x){
    nls.out <- nls.lm(par = list(rheobase = min(df[,x],na.rm=T)), fn = expresidFun, 
                      control = nls.lm.control(maxiter = 1024, maxfev = 1024), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
    print(nls.out)
    unlist(nls.out$par[1])
  }
  
  result <- Vectorize(fitDecay)(1:ncol(df))
  print(result)
  
  output <- t(matrix(c(colnames(df),result,rep(foldername,ncol(df))),nrow=3,byrow=T))
  colnames(output) <- c("Ind","rheobase","Genotype")
  all_fits <- rbind(all_fits,output)
}

all_fits <- as.data.frame(all_fits[-1,])

#write.csv(all_fits,"OutFiles/Rheobase/exp_decay_rheo_method.csv",row.names=F)


all_fits <- as.data.frame(all_fits)
comparison <- kruskal.test(all_fits$rheobase ~ all_fits$Genotype)
print(comparison)

# dunn <- c("exp_decay_t",dunn.test::dunn.test(all_fits$t, g = all_fits$Genotype, 
#                                        kw = T, method = "bonferroni",
#                                        table = T, list = F))

dunn.test::dunn.test(as.numeric(all_fits$rheobase), g = all_fits$Genotype, 
                     kw = T, method = "bonferroni",
                     table = T, list = F)
