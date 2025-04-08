#Libraries ----- 
library(ggplot2)
library(tidyverse)

#Aggregate C4P Correlation Results ----
c4p <- read.csv("OutFiles/C4P/test/Couchii_C4P_Metrics.csv")
head(c4p)

#set up matrix to store results
c4p_corr <- matrix(nrow=22,ncol=8)
colnames(c4p_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (col in 7:17) {
  pearson_corr <- cor.test(c4p$MAMU, c4p[,col], method = "pearson")
  #print(colnames(c4p[col])) here to show when the loop breaks

  c4p_corr[index,] <- c(colnames(c4p[col]),"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate)
  
  index = index + 1
  
  kendall_corr <- cor.test(c4p$MAMU, c4p[,col], method = "kendall")

  c4p_corr[index,] <- c(colnames(c4p[col]),"Kendall",kendall_corr$statistic,
                            NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
  
  index = index + 1
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_corr)) {
  if (c4p_corr[i,5] < 0.05) {
    print(c4p_corr[i,c(1,2,3,5)])
  }
}

#**I want to save this table as a file, with p<0.05 indicated in a column somehow---

#Splitting c4p data by pulse ----
c4p_pulse1 <- c4p[c4p$Pulse == 1,]
c4p_pulse2 <- c4p[c4p$Pulse == 2,]
c4p_pulse3 <- c4p[c4p$Pulse == 3,]
c4p_pulse4 <- c4p[c4p$Pulse == 4,]

#Comparing between pulses out of curiosity ----
#make long version of data:
c4p_long <- c4p %>% 
  select(Snake, Pulse:DiffFChgMaxToMin.ms.) %>% 
  pivot_longer(cols = BaseF.N.g.:DiffFChgMaxToMin.ms., 
               names_to = "variable", values_to = "value")
c4p_long$Pulse <- as.factor(c4p_long$Pulse)

pulse_boxes <- ggplot(c4p_long, aes(Pulse, value, fill = Pulse)) +
  geom_boxplot(outlier.shape = NA, na.rm = T) +
  facet_wrap (. ~ variable, scales = 'free', shrink = T) +
  xlab('') +
  ylab('')
pulse_boxes
#all this tells me is that at least visually all the pulses are the same, so I won't waste time comparing 
# differences between pulses between MAMUs or IC50s

#MAMU-Split pulse correlation results ----
c4p_corr_split <- matrix(nrow=88,ncol=9)
colnames(c4p_corr_split) <- c("Pulse","Param","CorrTest","Statistic", "df","p-value",
                        "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (pulse in 1:4) {
  df <- c4p[c4p$Pulse == pulse,]
  for (col in 7:17) {
    pearson_corr <- cor.test(df$MAMU, df[,col], method = "pearson")
    #print(colnames(c4p[col])) here to show when the loop breaks
    
    c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Pearson",pearson_corr$statistic,
                          pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                          pearson_corr$conf.int[2],pearson_corr$estimate)
    
    index = index + 1
    
    kendall_corr <- cor.test(df$MAMU, df[,col], method = "kendall")
    
    c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Kendall",kendall_corr$statistic,
                                NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
    
    index = index + 1
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_corr_split)) {
  if (c4p_corr_split[i,6] < 0.05) {
    print(c4p_corr_split[i,c(1,2,3,4,6)])
  }
}
#MAMU-split pulse linear regressions -----
par(mfrow=c(2,2))

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way

for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p[c4p$Pulse == pulse,]
    df[complete.cases(df[,18]),]
    plot(df$MAMU,df[,col],main = paste0("Pulse", pulse, colnames(df)[col]), 
         xlab = "TTX Resistance (MAMU)",
         ylab = colnames(df)[col])
    model <- lm(df[,col] ~ df$MAMU)
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
}

#**want to store all the linear regression results somehow ----

#IC50-Split pulse correlation results ----
IC50c4p_corr_split <- matrix(nrow=88,ncol=9)
colnames(IC50c4p_corr_split) <- c("Pulse","Param","CorrTest","Statistic", "df","p-value",
                              "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (pulse in 1:4) {
  df <- c4p[c4p$Pulse == pulse,]
  df <- merge(df,IC50, by="Snake")
  df <- df[!is.na(df$IC50),]
  for (col in 7:17) {
    pearson_corr <- cor.test(df$IC50, df[,col], method = "pearson")
    #print(colnames(c4p[col])) here to show when the loop breaks
    
    IC50c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Pearson",pearson_corr$statistic,
                                pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                                pearson_corr$conf.int[2],pearson_corr$estimate)
    
    index = index + 1
    
    kendall_corr <- cor.test(df$IC50, df[,col], method = "kendall")
    
    IC50c4p_corr_split[index,] <- c(pulse, colnames(df[col]),"Kendall",kendall_corr$statistic,
                                NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
    
    index = index + 1
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(IC50c4p_corr_split)) {
  if (IC50c4p_corr_split[i,6] < 0.05) {
    print(IC50c4p_corr_split[i,c(1,2,3,4,6)])
  }
}
#wait there are actually significant results here
#FMaxRate and FMinRate for all 4 pulses, Kendall only

#IC50-split pulse linear regressions -----
par(mfrow=c(2,2))

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way

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
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                    r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    abline(model, lwd=2, col="darkred")
    legend(x = "bottomright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
  }
}

#Repeat all c4p analyses without outliers ----
#(just ignoring this half at this point in my tidying and eventually will delete it 
# bc outliers will be removed at the end of data processing)

#remove outliers based on visual waveforms checks
c4p_sub <- c4p %>% 
  filter(!Snake %in% c("CRF3066", "CRF2677", "CRF2671", "CRF2669"))

#set up matrix to store results
c4p_sub_corr <- matrix(nrow=22,ncol=8)
colnames(c4p_sub_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                        "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (col in 7:17) {
  pearson_corr <- cor.test(c4p_sub$MAMU, c4p_sub[,col], method = "pearson")
  #print(colnames(c4p_sub[col])) here to show when the loop breaks
  
  c4p_sub_corr[index,] <- c(colnames(c4p_sub[col]),"Pearson",pearson_corr$statistic,
                        pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                        pearson_corr$conf.int[2],pearson_corr$estimate)
  
  index = index + 1
  
  kendall_corr <- cor.test(c4p_sub$MAMU, c4p_sub[,col], method = "kendall")
  
  c4p_sub_corr[index,] <- c(colnames(c4p_sub[col]),"Kendall",kendall_corr$statistic,
                        NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
  
  index = index + 1
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_sub_corr)) {
  if (c4p_sub_corr[i,5] < 0.05) {
    print(c4p_sub_corr[i,c(1,2,3,5)])
  }
}

#Splitting c4p_sub data by pulse 
c4p_sub_pulse1 <- c4p_sub[c4p_sub$Pulse == 1,]
c4p_sub_pulse2 <- c4p_sub[c4p_sub$Pulse == 2,]
c4p_sub_pulse3 <- c4p_sub[c4p_sub$Pulse == 3,]
c4p_sub_pulse4 <- c4p_sub[c4p_sub$Pulse == 4,]

#Comparing between pulses out of curiosity 
#make long version of data:
c4p_sub_long <- c4p_sub %>% 
  select(Snake, Pulse:DiffFChgMaxToMin.ms.) %>% 
  pivot_longer(cols = BaseF.N.g.:DiffFChgMaxToMin.ms., 
               names_to = "variable", values_to = "value")
c4p_sub_long$Pulse <- as.factor(c4p_sub_long$Pulse)

pulse_boxes <- ggplot(c4p_sub_long, aes(Pulse, value, fill = Pulse)) +
  geom_boxplot(outlier.shape = NA, na.rm = T) +
  facet_wrap (. ~ variable, scales = 'free', shrink = T) +
  xlab('') +
  ylab('')
pulse_boxes
#all this tells me is that at least visually all the pulses are the same, so I won't waste time comparing 
# differences between pulses between MAMUs or IC50s

#MAMU-Split pulse correlation results 
c4p_sub_corr_split <- matrix(nrow=88,ncol=9)
colnames(c4p_sub_corr_split) <- c("Pulse","Param","CorrTest","Statistic", "df","p-value",
                              "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (pulse in 1:4) {
  df <- c4p_sub[c4p_sub$Pulse == pulse,]
  for (col in 7:17) {
    pearson_corr <- cor.test(df$MAMU, df[,col], method = "pearson")
    #print(colnames(c4p_sub[col])) here to show when the loop breaks
    
    c4p_sub_corr_split[index,] <- c(pulse, colnames(df[col]),"Pearson",pearson_corr$statistic,
                                pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                                pearson_corr$conf.int[2],pearson_corr$estimate)
    
    index = index + 1
    
    kendall_corr <- cor.test(df$MAMU, df[,col], method = "kendall")
    
    c4p_sub_corr_split[index,] <- c(pulse, colnames(df[col]),"Kendall",kendall_corr$statistic,
                                NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
    
    index = index + 1
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_sub_corr_split)) {
  if (c4p_sub_corr_split[i,6] < 0.05) {
    print(c4p_sub_corr_split[i,c(1,2,3,4,6)])
  }
}
#MAMU-split pulse linear regressions 
par(mfrow=c(2,2))

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way

for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p_sub[c4p_sub$Pulse == pulse,]
    df[complete.cases(df[,18]),]
    plot(df$MAMU,df[,col],main = paste0("Pulse", pulse, colnames(df)[col]), 
         xlab = "TTX Resistance (MAMU)",
         ylab = colnames(df)[col])
    model <- lm(df[,col] ~ df$MAMU)
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
}

#IC50-Split pulse correlation results 
IC50c4p_sub_corr_split <- matrix(nrow=88,ncol=9)
colnames(IC50c4p_sub_corr_split) <- c("Pulse","Param","CorrTest","Statistic", "df","p-value",
                                  "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (pulse in 1:4) {
  df <- c4p_sub[c4p_sub$Pulse == pulse,]
  df <- merge(df,IC50, by="Snake")
  df <- df[!is.na(df$IC50),]
  for (col in 7:17) {
    pearson_corr <- cor.test(df$IC50, df[,col], method = "pearson")
    #print(colnames(c4p_sub[col])) here to show when the loop breaks
    
    IC50c4p_sub_corr_split[index,] <- c(pulse, colnames(df[col]),"Pearson",pearson_corr$statistic,
                                    pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                                    pearson_corr$conf.int[2],pearson_corr$estimate)
    
    index = index + 1
    
    kendall_corr <- cor.test(df$IC50, df[,col], method = "kendall")
    
    IC50c4p_sub_corr_split[index,] <- c(pulse, colnames(df[col]),"Kendall",kendall_corr$statistic,
                                    NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
    
    index = index + 1
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(IC50c4p_sub_corr_split)) {
  if (IC50c4p_sub_corr_split[i,6] < 0.05) {
    print(IC50c4p_sub_corr_split[i,c(1,2,3,4,6)])
  }
}
#wait there are actually significant results here, but different from the old sig results

#IC50-split pulse linear regressions 
par(mfrow=c(2,2))

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way

for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p_sub[c4p_sub$Pulse == pulse,]
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
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                    r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    abline(model, lwd=2, col="darkred")
    legend(x = "bottomright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
  }
}
