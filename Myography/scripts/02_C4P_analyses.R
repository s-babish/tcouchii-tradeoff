#C4P - Data analysis
# Computes Pearson's correlation coefficient and Kendall's Tau for MAMU and IC50
# for all the transient contraction summary stats computed by script 01. Also 
# does basic linear regressions of the summary stats on MAMU and IC50. All
# correlation coefficients and linear regression parameters (and r-squared and RMSE)
# are written to .csvs

#Libraries ----- 
library(tidyverse)

#Aggregate C4P Correlation Results ----
c4p <- read.csv("OutFiles/C4P/test/Couchii_C4P_Metrics.csv")
head(c4p)

#commented out because even if the pulses aren't distinguishable I'm not sure if
# this is valid (and I feel like I should do one or the other but not both)

# #set up matrix to store results
# c4p_corr <- matrix(nrow=22,ncol=9)
# colnames(c4p_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
#                             "Conf_int_LB","Conf_int_HB","CorEst","Sig")
# 
# index = 1
# for (col in 7:17) {
#   pearson_corr <- cor.test(c4p$MAMU, c4p[,col], method = "pearson")
# 
#   c4p_corr[index,] <- c(colnames(c4p[col]),"Pearson",pearson_corr$statistic,
#                             pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
#                             pearson_corr$conf.int[2],pearson_corr$estimate,"N")
#   
#   index = index + 1
#   
#   kendall_corr <- cor.test(c4p$MAMU, c4p[,col], method = "kendall")
# 
#   c4p_corr[index,] <- c(colnames(c4p[col]),"Kendall",kendall_corr$statistic,
#                             NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate,"N")
#   
#   index = index + 1
# }
# 
# #make it print out the ones with p<0.05
# for (i in 1:nrow(c4p_corr)) {
#   if (c4p_corr[i,5] < 0.05) {
#     print(c4p_corr[i,c(1,2,3,5)])
#     c4p_corr[i,9] <- "Y"
#   }
# }
# 
# write.csv(c4p_corr,"OutFiles/C4P/test/Couchii_C4P_MAMU_corr.csv")

#MAMU-Split pulse correlation results ----
c4p_corr_split <- matrix(nrow=88,ncol=10)
colnames(c4p_corr_split) <- c("Pulse","Param","CorrTest","Statistic", "df","p-value",
                        "Conf_int_LB","Conf_int_HB","CorEst", "Sig")

index = 1
for (pulse in 1:4) {
  df <- c4p[c4p$Pulse == pulse,]
  for (col in 7:17) {
    pearson_corr <- cor.test(df$MAMU, df[,col], method = "pearson")
    
    c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Pearson",pearson_corr$statistic,
                          pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                          pearson_corr$conf.int[2],pearson_corr$estimate,"N")
    
    index = index + 1
    
    kendall_corr <- cor.test(df$MAMU, df[,col], method = "kendall")
    
    c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Kendall",kendall_corr$statistic,
                                NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate,"N")
    
    index = index + 1
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_corr_split)) {
  if (c4p_corr_split[i,6] < 0.05) {
    print(c4p_corr_split[i,c(1,2,3,4,6)])
    c4p_corr_split[i,10] <- "Y"
  }
}

write.csv(c4p_corr,"OutFiles/C4P/test/Couchii_C4P_MAMU_corr_split.csv")

#MAMU-split pulse linear regressions -----
par(mfrow=c(2,2))

#make a storage matrix to keep the coefficients and RMSE/R^2 in
c4p_MAMU_reg <- matrix(nrow = 44,ncol = 6)
colnames(c4p_MAMU_reg) <- c("Metric","Pulse","RMSE","R2","B0","B1")

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way
i=1
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p[c4p$Pulse == pulse,]
    df[complete.cases(df[,18]),]
    plot(df$MAMU,df[,col],main = paste0("Pulse", pulse, colnames(df)[col]), 
         xlab = "TTX Resistance (MAMU)",
         ylab = colnames(df)[col])
    
    #do the actual linear regression
    model <- lm(df[,col] ~ df$MAMU)
    rmse <- round(sqrt(mean(resid(model)^2)), 2)
    coefs <- coef(model)
    b0 <- round(coefs[1], 2)
    b1 <- round(coefs[2],2)
    r2 <- round(summary(model)$r.squared, 2)
    
    #save all the values we just calculated
    c4p_MAMU_reg[i,] <- c(colnames(df)[col],pulse,rmse,r2,b0,b1)
    
    #now add them to the plot
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                    r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    abline(model, lwd=2, col="darkred")
    legend(x = "bottomright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
    i=i+1
  }
}

#save the storage matrix
write.csv(c4p_MAMU_reg,"OutFiles/C4P/test/Couchii_C4P_MAMU_lm.csv")

#IC50-Split pulse correlation results ----
IC50 <- read.csv("data_raw/IC50/IC50.csv")
IC50c4p_corr_split <- matrix(nrow=88,ncol=10)
colnames(IC50c4p_corr_split) <- c("Pulse","Param","CorrTest","Statistic", "df","p-value",
                              "Conf_int_LB","Conf_int_HB","CorEst", "Sig")

index = 1
for (pulse in 1:4) {
  df <- c4p[c4p$Pulse == pulse,]
  df <- merge(df,IC50, by="Snake")
  df <- df[!is.na(df$IC50),]
  for (col in 7:17) {
    pearson_corr <- cor.test(df$IC50, df[,col], method = "pearson")
    
    IC50c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Pearson",pearson_corr$statistic,
                                pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                                pearson_corr$conf.int[2],pearson_corr$estimate, "N")
    
    index = index + 1
    
    kendall_corr <- cor.test(df$IC50, df[,col], method = "kendall")
    
    IC50c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Kendall",kendall_corr$statistic,
                                NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate, "N")
    
    index = index + 1
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(IC50c4p_corr_split)) {
  if (IC50c4p_corr_split[i,6] < 0.05) {
    print(IC50c4p_corr_split[i,c(1,2,3,4,6)])
    IC50c4p_corr_split[i,10] <- "Y"
  }
}

write.csv(IC50c4p_corr_split,"OutFiles/C4P/test/Couchii_C4P_IC50_corr_split.csv")

#IC50-split pulse linear regressions -----
par(mfrow=c(2,2))

#storage matrix again
c4p_IC50_reg <- matrix(nrow = 44,ncol = 6)
colnames(c4p_IC50_reg) <- c("Metric","Pulse","RMSE","R2","B0","B1")

j=1
#same process as for the MAMUs
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p[c4p$Pulse == pulse,]
    df <- merge(df,IC50, by="Snake")
    df <- df[!is.na(df$IC50),]
    plot(df$IC50,df[,col],main = paste0("Pulse", pulse, colnames(df)[col]), 
         xlab = "TTX Resistance (IC50)",
         ylab = colnames(df)[col])
    
    model <- lm(df[,col] ~ df$IC50)
    rmse <- round(sqrt(mean(resid(model)^2)), 2)
    coefs <- coef(model)
    b0 <- round(coefs[1], 2)
    b1 <- round(coefs[2],2)
    r2 <- round(summary(model)$r.squared, 2)
    c4p_IC50_reg[j,] <- c(colnames(df)[col],pulse,rmse,r2,b0,b1)
    
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                    r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    abline(model, lwd=2, col="darkred")
    legend(x = "bottomright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
    
    j=j+1
  }
}

#storage file yet again
write.csv(c4p_IC50_reg,"OutFiles/C4P/test/Couchii_C4P_IC50_lm.csv")
