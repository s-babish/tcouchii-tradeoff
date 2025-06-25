#all old split-pulse analyses got moved into here
#MAMU-Split pulse correlation results ----
c4p_corr_split <- matrix(nrow=44,ncol=10)
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
  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_corr_split)) {
  if (c4p_corr_split[i,6] < 0.05) {
    print(c4p_corr_split[i,c(1,2,3,4,6)])
    c4p_corr_split[i,10] <- "Y"
  }
}

#write.csv(c4p_corr_split,"OutFiles/C4P/Couchii_C4P_MAMU_corr_split.csv")

#MAMU-split pulse linear regressions -----
par(mfrow=c(2,2))

#Check for normality in data subset and transform as needed 
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p[c4p$Pulse == pulse,]
    df[complete.cases(df[,18]),]
    print(pulse)
    print(colnames(df)[col])
    print(shapiro.test(df[,col]))
  }
}

#FMinRateOfCHg, FMaxRateOfChg, ContrAmpl, and BaseF not normal
# cols 7, 8, 13, 14 need log transformed

c4p_log <- c4p
logcols <- c(7,8,13,14)

for (col in logcols) {

  c4p_log[,col] <- log(abs(c4p[,col]))
  colnames(c4p_log)[col] <- paste0("Log",colnames(c4p)[col])
}

#re-check for normality
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p_log[c4p_log$Pulse == pulse,]
    df[complete.cases(df[,18]),]
    print(pulse)
    print(colnames(df)[col])
    print(shapiro.test(df[,col]))
  }
}

#write function to extract p value
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#make a storage matrix to keep the coefficients and RMSE/R^2 in
c4p_MAMU_reg <- matrix(nrow = 44,ncol = 7)
colnames(c4p_MAMU_reg) <- c("Metric","Pulse","RMSE","R2","B0","B1", "p")

#this one will do it all of the same value for each pulse simultaneously
#i know it's less efficient code-wise but it's easier to compare this way
i=1
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p_log[c4p_log$Pulse == pulse,]
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
    p <- overall_p(model)
    #save all the values we just calculated
    c4p_MAMU_reg[i,] <- c(colnames(df)[col],pulse,rmse,r2,b0,b1,p)

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
#write.csv(c4p_MAMU_reg,"OutFiles/C4P/Couchii_C4P_MAMU_lm.csv")

#IC50-Split pulse correlation results ----
IC50 <- read.csv("data_raw/IC50/IC50.csv")
IC50c4p_corr_split <- matrix(nrow=44,ncol=10)
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

  }
}

#make it print out the ones with p<0.05
for (i in 1:nrow(IC50c4p_corr_split)) {
  if (IC50c4p_corr_split[i,6] < 0.05) {
    print(IC50c4p_corr_split[i,c(1,2,3,4,6)])
    IC50c4p_corr_split[i,10] <- "Y"
  }
}

#write.csv(IC50c4p_corr_split,"OutFiles/C4P/Couchii_C4P_IC50_corr_split.csv")

#IC50-split pulse linear regressions -----
par(mfrow=c(2,2))

#Check for normality in data subset and transform as needed ----
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p[c4p$Pulse == pulse,]
    df <- merge(df,IC50, by="Snake")
    df <- df[!is.na(df$IC50),]

    print(pulse)
    print(colnames(df)[col])
    print(shapiro.test(df[,col]))
  }
}

#BaseF not normal
# col 7 needs log transformed
c4p_log[,7] <- log(c4p[,7])
colnames(c4p_log)[7] <- paste0("Log",colnames(c4p)[7])

#re-check for normality
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p_log[c4p_log$Pulse == pulse,]
    df <- merge(df,IC50, by="Snake")
    df <- df[!is.na(df$IC50),]

    print(pulse)
    print(colnames(df)[col])
    print(shapiro.test(df[,col]))
  }
}

#storage matrix again
c4p_IC50_reg <- matrix(nrow = 44,ncol = 7)
colnames(c4p_IC50_reg) <- c("Metric","Pulse","RMSE","R2","B0","B1","p")

j=1
#same process as for the MAMUs
for (col in 7:17) {
  for (pulse in 1:4) {
    df <- c4p_log[c4p_log$Pulse == pulse,]
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
    p <- overall_p(model)

    c4p_IC50_reg[j,] <- c(colnames(df)[col],pulse,rmse,r2,b0,b1,p)

    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~
                    r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    abline(model, lwd=2, col="darkred")
    legend(x = "bottomright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))
    #plot(model)
    #par(mfrow=c(1,1))
    j=j+1
  }
}

#storage file yet again
#write.csv(c4p_IC50_reg,"OutFiles/C4P/Couchii_C4P_IC50_lm.csv")
