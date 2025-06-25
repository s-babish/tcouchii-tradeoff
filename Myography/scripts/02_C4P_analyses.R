#C4P - Data analysis
# Computes Pearson's correlation coefficient and Kendall's Tau for MAMU and IC50
# for all the transient contraction summary stats computed by script 01. Also 
# does basic linear regressions of the summary stats on MAMU and IC50. All
# correlation coefficients and linear regression parameters (and r-squared and RMSE)
# are written to .csvs

#Libraries ----- 
library(tidyverse)
library(MASS)

#Load data ----
c4p <- read.csv("OutFiles/C4P/Couchii_C4P_Metrics.csv")
IC50 <- read.csv("data_raw/IC50/IC50.csv") #for later
#remove outliers ID'd at end of script 1
c4p <- c4p %>% 
  filter(!Snake %in% c("CRF3066", "CRF2677", "CRF2671", "CRF2669"))

#determine which data is normal and transform whatever isn't ----
nonnormal <- numeric() #storage vector for indices of non-normal columns
for (col in 7:17) {
  normality_test <- shapiro.test(c4p[,col])
  if(normality_test$p.value < 0.05) {
    print((colnames(c4p)[col])) #print names of columns not normally distributed
    nonnormal <- c(nonnormal,col)
  }
}
#not literally all of them not being normal
log_c4p <- c4p
log_c4p[,c(7:17)] <- log(abs(log_c4p[,c(7:17)])) #have to take abs value by FMinRate is negative

nonnormal <- numeric() #storage vector for indices of non-normal columns
for (col in 7:17) {
  normality_test <- shapiro.test(log_c4p[,col])
  if(normality_test$p.value < 0.05) {
    print((colnames(log_c4p)[col])) #print names of columns not normally distributed
    nonnormal <- c(nonnormal,col)
  }
}

#option 1:: kruskal-wallis test (as non-parametric anova equivalent) ----
for (col in 7:17) {
  print(colnames(c4p)[col])
  test <- kruskal.test(c4p[,col] ~ c4p$Pulse)
  print(test)
  # all pulses are comprehensively identical, move on
}

#option 2: box-cox transformation and then anova ----
bc_lambdas <- numeric()
bc_c4p <- c4p #set up to transform cols
for (col in c(7,8,9,10,11,12,13,15,16,17)) { #can't use for negative values
  x <- c4p[,col]
  bc <- boxcox(lm(x~1), lambda = seq(-5, 5, 1/10)) #i'm using full range of allowed values bc why not
  lambda <- bc$x[which.max(bc$y)] #pull out optimal lambda value
  bc_lambdas <- c(bc_lambdas, lambda) #save new optimal lambda
  
  bc_c4p[,col] <- (x ^ lambda - 1) / lambda
  colnames(bc_c4p)[col] <- paste0(colnames(c4p)[col],"boxcox",lambda)
}
#not going to actually use this now beyond the anova, but will if i come up with something applicable
# (maybe for linear regressions? not correlations though)

for(col in 7:17) {
  print(colnames(bc_c4p)[col])
  anova <- aov(bc_c4p[,col] ~ bc_c4p$Pulse)
  print(summary(anova))
}
#this also confirmed pulses were different, hooray i can move on now

#(Aggregate) C4P MAMU Correlation Results ----
#set up matrix to store results
c4p_corr <- matrix(nrow=11,ncol=9)
colnames(c4p_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                        "Conf_int_LB","Conf_int_HB","CorEst","Sig")

index = 1
for (col in 7:17) {
  pearson_corr <- cor.test(c4p$MAMU, c4p[,col], method = "pearson")
  
  c4p_corr[index,] <- c(colnames(c4p)[col],"Pearson",pearson_corr$statistic,
                        pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                        pearson_corr$conf.int[2],pearson_corr$estimate,"N")
  
  index = index + 1
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_corr)) {
  if (c4p_corr[i,5] < 0.05) {
    print(c4p_corr[i,c(1,2,3,5)])
    c4p_corr[i,9] <- "Y"
  }
}

write.csv(c4p_corr,"OutFiles/C4P/Couchii_C4P_MAMU_corr.csv")

#C4P MAMU regressions ----
#write function to extract p value
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#make a storage matrix to keep the coefficients and RMSE/R^2 in
c4p_MAMU_reg <- matrix(nrow = 11,ncol = 6)
colnames(c4p_MAMU_reg) <- c("Metric","RMSE","R2","B0","B1", "p")

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way

df<-bc_c4p #lazy and don't feel like changing all the df's below
i=1
for (col in 7:17) {
  df[complete.cases(df[,18]),]
  plot(df$MAMU,df[,col],main = paste0(colnames(df)[col]),
       xlab = "TTX Resistance (MAMU)",
       ylab = colnames(df)[col])
  
  #do the actual linear regression
  model <- lm(df[,col] ~ df$MAMU)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  p <- overall_p(model)
  #save all the values we just calculated
  c4p_MAMU_reg[i,] <- c(colnames(df)[col],rmse,r2,b0,b1,p)
  
  #now add them to the plot
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
  i=i+1
}


#save the storage matrix
write.csv(c4p_MAMU_reg,"OutFiles/C4P/Couchii_C4P_MAMU_lm.csv")

#C4P IC50 correlations ----
#make storage matrix
IC50c4p_corr <- matrix(nrow=11,ncol=9)
colnames(IC50c4p_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst", "Sig")

df <- merge(c4p,IC50, by="Snake")
df <- df[!is.na(df$IC50),]

index = 1
for (col in 7:17) {
  pearson_corr <- cor.test(df$IC50, df[,col], method = "pearson")
  
  IC50c4p_corr[index,] <- c(colnames(df)[col],"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate, "N")
  index = index + 1
}


#make it print out the ones with p<0.05
for (i in 1:nrow(IC50c4p_corr)) {
  if (IC50c4p_corr[i,5] < 0.05) {
    print(IC50c4p_corr[i,c(1,5,8)])
    IC50c4p_corr[i,9] <- "Y"
  }
}

write.csv(IC50c4p_corr_split,"OutFiles/C4P/Couchii_C4P_IC50_corr.csv")

#C4P IC50 regressions -----
df <- merge(bc_c4p,IC50, by = "Snake")
df <- df[!is.na(df$IC50),]

#storage matrix again
c4p_IC50_reg <- matrix(nrow = 11,ncol = 6)
colnames(c4p_IC50_reg) <- c("Metric","RMSE","R2","B0","B1","p")

index=1
#same process as for the MAMUs
for (col in 7:17) {
    plot(df$IC50,df[,col],main = paste0(colnames(df)[col]),
         xlab = "TTX Resistance (IC50)",
         ylab = colnames(df)[col])
    
    model <- lm(df[,col] ~ df$IC50)
    rmse <- round(sqrt(mean(resid(model)^2)), 2)
    coefs <- coef(model)
    b0 <- round(coefs[1], 2)
    b1 <- round(coefs[2],2)
    r2 <- round(summary(model)$r.squared, 2)
    p <- overall_p(model)
    
    c4p_IC50_reg[index,] <- c(colnames(df)[col],rmse,r2,b0,b1,p)
    
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~
                    r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    abline(model, lwd=2, col="darkred")
    legend(x = "bottomright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
    #plot(model)
    #par(mfrow=c(1,1))
    index=index+1
  
}

#storage file yet again
write.csv(c4p_IC50_reg,"OutFiles/C4P/Couchii_C4P_IC50_lm.csv")

#rerun the transformed regressions with just 1 pulse involved to see the effect 
# of removing the repeats ----
df<-df[df$Pulse == 1,]
c4p_IC50_reg <- matrix(nrow = 11,ncol = 6)
colnames(c4p_IC50_reg) <- c("Metric","RMSE","R2","B0","B1","p")

index=1
#same process as for the MAMUs
for (col in 7:17) {
  plot(df$IC50,df[,col],main = paste0(colnames(df)[col]),
       xlab = "TTX Resistance (IC50)",
       ylab = colnames(df)[col])
  
  model <- lm(df[,col] ~ df$IC50)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  p <- overall_p(model)
  
  c4p_IC50_reg[index,] <- c(colnames(df)[col],rmse,r2,b0,b1,p)
  
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
  #plot(model)
  #par(mfrow=c(1,1))
  index=index+1
}

#that did remove all the correlations