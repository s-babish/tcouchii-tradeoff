setwd("C:/Users/sdbab/OneDrive - University of Nevada, Reno/UNR/tcouchii/Myography")

#Libraries ----- 
library(ggplot2)
library(tidyverse)


#Aggregate C4P Correlation Results ----
c4p <- read.csv("Tcouchii_side_hustle/Couchii_C4P/output/Couchii_C4P.csv")
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

#MAMU-Tetanus correlation results ----
tetanus <- read.csv("Tcouchii_side_hustle/Couchii_tetanus/output/Couchii_tetanus.csv")
#remove the duplicate rows that are in here for some reason (note that this method only works bc
# there's only one muscle from each snake, I checked):
tetanus <- tetanus[!duplicated(tetanus$Snake), ]
head(tetanus)

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



#Rheobase correlation results (kind of) -----
#doesn't work rn
pulse_l <- c(50,100,200,500,1000,10000,50000)

for (i in 1:length(pulse_l)) {
  infile <- toString(paste(pulse_l[i],"_TPA-sigmoidal-rpt.csv", sep=""))
  df <- read.csv(infile)
  par(mfrow=c(1,3))
  for (col in 3:5) {
    print(pulse_l, colnames(df[col]))
    
    pearson_corr <- cor.test(df$MAMU,df[,col], method = "pearson")
    print(pearson_corr)
    
    kendall_corr <- cor.test(df$MAMU,df[,col], method = "kendall")
    print(kendall_corr)
  }
}


