#want to implement all the changes/new data saves I made in the C4p analyses here
#(i can see a lot of them ending up as supplementals)

#Libraries ----- 
library(ggplot2)
library(tidyverse)

#MAMU-Tetanus correlation results ----
tetanus <- read.csv("OutFiles/Tetanus/test/Couchii_Tetanus_Metrics.csv")
#remove the duplicate rows that are in here for some reason (note that this method only works bc
# there's only one muscle from each snake, I checked):
tetanus <- tetanus[!duplicated(tetanus$Snake), ] 

#set up matrix to store results
tetanus_corr <- matrix(nrow=20,ncol=9)
colnames(tetanus_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst","Sig")
index = 1
for (col in 4:13) {
  pearson_corr <- cor.test(tetanus$MAMU,tetanus[,col], method = "pearson")
  tetanus_corr[index,] <- c(colnames(tetanus[col]),"Pearson",pearson_corr$statistic,
                        pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                        pearson_corr$conf.int[2],pearson_corr$estimate,"N")
                        
  index = index + 1
  
  kendall_corr <- cor.test(tetanus$MAMU,tetanus[,col], method = "kendall")
  tetanus_corr[index,] <- c(colnames(tetanus[col]),"Kendall",kendall_corr$statistic,
                           NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate,"N")
                             
  index = index + 1
}

#make it print the ones with p<0.05
for (i in 1:nrow(tetanus_corr)) {
  if (tetanus_corr[i,5] < 0.05) {
    print(tetanus_corr[i,c(1,2,3,5)])
    tetanus_crr[i,9] <- "Y"
  }
}

write.csv(tetanus_corr, "OutFiles/Tetanus/test/Couchii_Tetanus_MAMU_corr.csv")

#MAMU-tetanus linear regressions -----

tetanus_MAMU_reg <- matrix(nrow = 20, ncol = 5)
colnames(tetanus_MAMU_reg) <- c("Metric","RMSE","R2","B0","B1")

i=1
for (col in 4:13) {
  plot(tetanus$MAMU,tetanus[,col],main = colnames(tetanus)[col], xlab = "TTX Resistance (MAMU)",
       ylab = colnames(tetanus)[col])
  
  model <- lm(tetanus[,col] ~ tetanus$MAMU)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  tetanus_MAMU_reg[i,] <- c(colnames(tetanus)[col],rmse,r2,b0,b1)
  
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
  i=i+1
}

write.csv(tetanus_MAMU_reg,"OutFiles/Tetanus/test/Couchii_Tetanus_MAMU_lm.csv")

#Correlation with IC50----
#read in IC50 and merge dataframes
#IC50tet <- tetanus
IC50 <- read.csv("data_raw/IC50/IC50.csv")
IC50tet <- merge(tetanus,IC50, by="Snake")
IC50tet <- IC50tet[!is.na(IC50tet$IC50),]

summary(IC50tet)
#only 15 obs in this comparison

#correlation metrics
IC50tet_corr <- matrix(nrow=20,ncol=9)
colnames(IC50tet_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst","Sig")
index = 1
for (col in 4:13) {
  #print(colnames(IC50tet[col])) here to show when the loop breaks
  pearson_corr <- cor.test(IC50tet$MAMU,IC50tet[,col], method = "pearson")
  IC50tet_corr[index,] <- c(colnames(IC50tet[col]),"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate,"N")
  
  index = index + 1
  
  kendall_corr <- cor.test(IC50tet$MAMU,IC50tet[,col], method = "kendall")
  IC50tet_corr[index,] <- c(colnames(IC50tet[col]),"Kendall",kendall_corr$statistic,
                            NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate,"N")
  
  index = index + 1
}

#make it print the ones with p<0.05
for (i in 1:nrow(IC50tet_corr)) {
  if (IC50tet_corr[i,5] < 0.05) {
    print(IC50tet_corr[i,c(1,2,3,5)])
    IC50tet_corr[i,9] <- "Y"
  }
}

write.csv(IC50tet_corr,"OutFiles/Tetanus/test/Couchii_Tetanus_IC50_corr.csv")

#IC50-tetanus linear regressions ----

#storage matrix again
tetanus_IC50_reg <- matrix(nrow = 20,ncol = 5)
colnames(tetanus_IC50_reg) <- c("Metric","RMSE","R2","B0","B1")

j=1
for (col in 4:13) {
  plot(IC50tet$IC50,IC50tet[,col],main = colnames(IC50tet)[col], xlab = " Muscle TTX Resistance/IC50 (MAMU)",
       ylab = colnames(IC50tet)[col])
  
  model <- lm(IC50tet[,col] ~ IC50tet$IC50)
  rmse <- round(sqrt(mean(resid(model)^2)), 2)
  coefs <- coef(model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(model)$r.squared, 2)
  tetanus_IC50_reg[j,] <- c(colnames(IC50tet)[col],rmse,r2,b0,b1)
  
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                  r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  abline(model, lwd=2, col="darkred")
  legend(x = "bottomright", bty = "n",
         legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
  
  j = j+1
}

#write to outfile
write.csv(tetanus_IC50_reg,"OutFiles/Tetanus/test/Couchii_Tetanus_IC50_lm.csv")
