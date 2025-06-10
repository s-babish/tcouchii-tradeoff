#C4P - Data analysis
# Computes Pearson's correlation coefficient and Kendall's Tau for MAMU and IC50
# for all the transient contraction summary stats computed by script 01. Also 
# does basic linear regressions of the summary stats on MAMU and IC50. All
# correlation coefficients and linear regression parameters (and r-squared and RMSE)
# are written to .csvs

#Libraries ----- 
library(tidyverse)

#Load data ----
c4p <- read.csv("OutFiles/C4P/Couchii_C4P_Metrics.csv")
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



#(Aggregate) C4P MAMU Correlation Results ----
#set up matrix to store results
c4p_corr <- matrix(nrow=22,ncol=9)
colnames(c4p_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                            "Conf_int_LB","Conf_int_HB","CorEst","Sig")

index = 1
for (col in 7:17) {
  pearson_corr <- cor.test(c4p$MAMU, c4p[,col], method = "pearson")

  c4p_corr[index,] <- c(colnames(c4p[col]),"Pearson",pearson_corr$statistic,
                            pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                            pearson_corr$conf.int[2],pearson_corr$estimate,"N")

  index = index + 1

  kendall_corr <- cor.test(c4p$MAMU, c4p[,col], method = "kendall")

  c4p_corr[index,] <- c(colnames(c4p[col]),"Kendall",kendall_corr$statistic,
                            NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate,"N")

  index = index + 1
}

#make it print out the ones with p<0.05
for (i in 1:nrow(c4p_corr)) {
  if (c4p_corr[i,5] < 0.05) {
    print(c4p_corr[i,c(1,2,3,5)])
    c4p_corr[i,9] <- "Y"
  }
}

write.csv(c4p_corr,"OutFiles/C4P/test/Couchii_C4P_MAMU_corr.csv")

