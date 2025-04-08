#want to implement all the changes/new data saves I made in the C4p analyses here
#(i can see a lot of them ending up as supplementals)

#Libraries ----- 
library(ggplot2)
library(tidyverse)

#MAMU-Tetanus correlation results ----
tetanus <- read.csv("OutFiles/Tetanus/test/Couchii_Tetanus_Metrics.csv")
#remove the duplicate rows that are in here for some reason (note that this method only works bc
# there's only one muscle from each snake, I checked):
tetanus_sub <- tetanus[!duplicated(tetanus$Snake), ] %>% 
  filter(!Snake %in% c("CRF3060","CRF3074","CRF3065","CRF3066","CRF2680","CRF2669","CRF2631","CRF2670")) 
#removing outliers based on waveform analysis currently in myography_plots.R (down to 24 obs)
head(tetanus_sub)

#set up matrix to store results
tetanus_corr <- matrix(nrow=20,ncol=8)
colnames(tetanus_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst")
index = 1
for (col in 4:13) {
  #print(colnames(tetanus[col])) here to show when the loop breaks
  pearson_corr <- cor.test(tetanus$MAMU,tetanus[,col], method = "pearson")
  tetanus_corr[index,] <- c(colnames(tetanus[col]),"Pearson",pearson_corr$statistic,
                        pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                        pearson_corr$conf.int[2],pearson_corr$estimate)
                        
  index = index + 1
  
  kendall_corr <- cor.test(tetanus$MAMU,tetanus[,col], method = "kendall")
  tetanus_corr[index,] <- c(colnames(tetanus[col]),"Kendall",kendall_corr$statistic,
                           NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
                             
  index = index + 1
}

#make it print the ones with p<0.05
for (i in 1:nrow(tetanus_corr)) {
  if (tetanus_corr[i,5] < 0.05) {
    print(tetanus_corr[i,c(1,2,3,5)])
  }
}

#MAMU-tetanus linear regressions -----
#linear regressions
for (col in 4:13) {
  plot(tetanus$MAMU,tetanus[,col],main = colnames(tetanus)[col], xlab = "TTX Resistance (MAMU)",
       ylab = colnames(tetanus)[col])
  model <- lm(tetanus[,col] ~ tetanus$MAMU)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
}

#Correlation with IC50tet (copied from IC50analyses.R)----
IC50tet <- tetanus
IC50tet <- merge(IC50tet,IC50, by="Snake")
IC50tet <- IC50tet[!is.na(IC50tet$IC50),]

summary(IC50tet)
#only 15 obs in this comparison

#correlation metrics
IC50tet_corr <- matrix(nrow=20,ncol=8)
colnames(IC50tet_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst")
index = 1
for (col in 4:13) {
  #print(colnames(IC50tet[col])) here to show when the loop breaks
  pearson_corr <- cor.test(IC50tet$MAMU,IC50tet[,col], method = "pearson")
  IC50tet_corr[index,] <- c(colnames(IC50tet[col]),"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate)
  
  index = index + 1
  
  kendall_corr <- cor.test(IC50tet$MAMU,IC50tet[,col], method = "kendall")
  IC50tet_corr[index,] <- c(colnames(IC50tet[col]),"Kendall",kendall_corr$statistic,
                            NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
  
  index = index + 1
}

#make it print the ones with p<0.05
for (i in 1:nrow(IC50tet_corr)) {
  if (IC50tet_corr[i,5] < 0.05) {
    print(IC50tet_corr[i,c(1,2,3,5)])
  }
}

#IC50-tetanus linear regressions ----
for (col in 4:13) {
  plot(IC50tet$IC50,IC50tet[,col],main = colnames(IC50tet)[col], xlab = " Muscle TTX Resistance/IC50 (MAMU)",
       ylab = colnames(IC50tet)[col])
  model <- lm(IC50tet[,col] ~ IC50tet$IC50)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
}

#Repeat all tetanus results with subset without outliers (separate for now to compare values) ----
#set up matrix to store results
tetanus_sub_corr <- matrix(nrow=20,ncol=8)
colnames(tetanus_sub_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst")
index = 1
for (col in 4:13) {
  #print(colnames(tetanus_sub[col])) here to show when the loop breaks
  pearson_corr <- cor.test(tetanus_sub$MAMU,tetanus_sub[,col], method = "pearson")
  tetanus_sub_corr[index,] <- c(colnames(tetanus_sub[col]),"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate)
  
  index = index + 1
  
  kendall_corr <- cor.test(tetanus_sub$MAMU,tetanus_sub[,col], method = "kendall")
  tetanus_sub_corr[index,] <- c(colnames(tetanus_sub[col]),"Kendall",kendall_corr$statistic,
                            NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
  
  index = index + 1
}

#make it print the ones with p<0.05
for (i in 1:nrow(tetanus_sub_corr)) {
  if (tetanus_sub_corr[i,5] < 0.05) {
    print(tetanus_sub_corr[i,c(1,2,3,5)])
  }
}

#MAMU-tetanus_sub linear regressions 
#linear regressions
for (col in 4:13) {
  plot(tetanus_sub$MAMU,tetanus_sub[,col],main = colnames(tetanus_sub)[col], xlab = "TTX Resistance (MAMU)",
       ylab = colnames(tetanus_sub)[col])
  model <- lm(tetanus_sub[,col] ~ tetanus_sub$MAMU)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
}

#Correlation with IC50tet (copied from IC50analyses.R)
IC50tet <- tetanus_sub
IC50tet <- merge(IC50tet,IC50, by="Snake")
IC50tet <- IC50tet[!is.na(IC50tet$IC50),]

summary(IC50tet)
#only 11 obs in this comparison

#correlation metrics
IC50tet_corr <- matrix(nrow=20,ncol=8)
colnames(IC50tet_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst")
index = 1
for (col in 4:13) {
  #print(colnames(IC50tet[col])) here to show when the loop breaks
  pearson_corr <- cor.test(IC50tet$MAMU,IC50tet[,col], method = "pearson")
  IC50tet_corr[index,] <- c(colnames(IC50tet[col]),"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate)
  
  index = index + 1
  
  kendall_corr <- cor.test(IC50tet$MAMU,IC50tet[,col], method = "kendall")
  IC50tet_corr[index,] <- c(colnames(IC50tet[col]),"Kendall",kendall_corr$statistic,
                            NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
  
  index = index + 1
}

#make it print the ones with p<0.05
for (i in 1:nrow(IC50tet_corr)) {
  if (IC50tet_corr[i,5] < 0.05) {
    print(IC50tet_corr[i,c(1,2,3,5)])
  }
}

#IC50-tetanus_sub linear regressions 
for (col in 4:13) {
  plot(IC50tet$IC50,IC50tet[,col],main = colnames(IC50tet)[col], xlab = " Muscle TTX Resistance/IC50 (MAMU)",
       ylab = colnames(IC50tet)[col])
  model <- lm(IC50tet[,col] ~ IC50tet$IC50)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
}
