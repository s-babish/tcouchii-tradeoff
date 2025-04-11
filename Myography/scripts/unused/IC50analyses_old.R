#IC50 expresses specifically the skeletal muscle resistance to TTX
#whole animal resistance and skeletal muscle resistance already known to be
# correlated (reimche et al. 2022) so i'm not sure what the point of this 
# analysis is

#used spearman's rank correlation which i might try to use for the other stuff
# just to be consistent (and I'll redo it here myself)

#data collected by dissecting out muscle and serially increasing TTX doses to
# see which reduced the peak transient force magnitude to 50%

#read in data
IC50 <- read.csv("IC50.csv")

head(IC50)

#remove individuals without IC50 data:

IC50 <- IC50[!(is.na(IC50$IC50)),]
dim(IC50)

#now have 25 observations; not the full 32/34 in the other datasets (c4p and 
# IC50tet and rheboase all vary)

#figure out which inds are unique to each dataset:

setdiff(c4p$Snake, tetanus$Snake)
#shared:
intersect(IC50$Collector,c4p$Snake)
#"CRF3051" "CRF3052" "CRF3056" "CRF3058" "CRF3059" "CRF3060" "CRF3061" "EJE164"  
#"CRF3070" "CRF3072" "CRF3074" "CRF3064" "CRF3065" "CRF3066" "EJE186"  "CRF3211"

#16 snakes overlap

#IC50 but not others:
setdiff(IC50$Collector, c4p$Snake)
#"CRF3062" "EJE163"  "EJE171"  "EJE166"  "CRF3293" "CRF3309" "CRF3221" "EJE182"  
#"CRF3217"

#9 snakes we only have IC50 for

#others but not IC50:
setdiff(c4p$Snake, IC50$Collector)
#"CRF2630" "CRF2631" "CRF2633" "CRF2669" "CRF2670" "CRF2671" "CRF2672" "CRF2673" 
#"CRF2674" "CRF2675"
#"CRF2676" "CRF2677" "CRF2678" "CRF2679" "CRF2680" "CRF2681" "CRF3055" "CRF3069"

#18 snakes in c4p but not both

#annoying because it makes relating the different results much harder

#Correlation Coefficients ----
#Spearman's correlation coefficient (p)
#assigns/assesses correlation based on rank
#good for monotonic relationships, so first will plot to confirm it's monotonic

plot(IC50$X50_MAMU, IC50$IC50)

#i guess this could be called monotonic? but it's also just vaguely linear
#i'd bet anything they used this one bc it gave a better number lol

spearman_IC50 <- cor.test(IC50$IC50,IC50$X50_MAMU, method = "spearman")
spearman_IC50
#also estimate with direct hypothesis (positive corr):
spearman_IC50_directional <- cor.test(IC50$IC50,IC50$X50_MAMU, method = "spearman",
                                     alternative = "greater")
spearman_IC50_directional
#marginally increases the p-value but it doesn't truly matter here

#Kendall's Tau
#after some reading I am also using this one because it can apparently be better
# for smaller sample sizes and if there are outliers

kendall_IC50 <- cor.test(IC50$IC50,IC50$X50_MAMU, method = "kendall")
kendall_IC50

#also estimate with direct hypothesis (positive corr):
kendall_IC50_directional <- cor.test(IC50$IC50,IC50$X50_MAMU, method = "kendall",
                                     alternative = "greater")
kendall_IC50_directional
#this doesn't change the test statistic but does increase the p-value slightly

#Linear Regression ----

par(mfrow=c(1,1))
IC50_lm <- lm(IC50$IC50 ~ IC50$X50_MAMU)
summary(IC50_lm)

rmse <- round(sqrt(mean(resid(IC50_lm)^2)), 2)
coefs <- coef(IC50_lm)
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)
r2 <- round(summary(IC50_lm)$r.squared, 2)
eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                r^2 == .(r2) * "," ~~ RMSE == .(rmse))
plot(IC50$X50_MAMU, IC50$IC50)
abline(IC50_lm, lwd=2, col="darkred")
legend(x = "bottomright", bty = "n",
       legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse)))

#wow a linear regression that actually makes sense
#not that the result makes sense logically but ya know


#Non-linear regressions (bc not normally distributed)

#Correlation with c4p ----
#add IC50 column to c4p dataset, and remove rows without both

IC50c4p <- merge(c4p,IC50, by="Snake")
IC50c4p <- IC50c4p[!is.na(IC50c4p$IC50),]

#linear regressions
for (col in 7:17) {
  plot(IC50c4p$IC50,IC50c4p[,col],main = colnames(IC50c4p)[col], xlab = "TTX Resistance (MAMU)",
       ylab = colnames(IC50c4p)[col])
  model <- lm(IC50c4p[,col] ~ IC50c4p$IC50)
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
#some slightly better r@ values: 0.22 for contrampl, 0.23 for fmaxrate ngs,
# 0.22 for fminrate ngs (but with bad RMSE)

#correlation analyses
IC50c4p_corr <- matrix(nrow=22,ncol=8)
colnames(IC50c4p_corr) <- c("Param","CorrTest","Statistic", "df","p-value",
                        "Conf_int_LB","Conf_int_HB","CorEst")

index = 1
for (col in 7:17) {
  pearson_corr <- cor.test(IC50c4p$IC50, IC50c4p[,col], method = "pearson")
  #print(colnames(IC50c4p[col])) here to show when the loop breaks
  
  IC50c4p_corr[index,] <- c(colnames(IC50c4p[col]),"Pearson",pearson_corr$statistic,
                        pearson_corr$parameter,pearson_corr$p.value,pearson_corr$conf.int[1],
                        pearson_corr$conf.int[2],pearson_corr$estimate)
  
  index = index + 1
  
  kendall_corr <- cor.test(IC50c4p$IC50, IC50c4p[,col], method = "kendall")
  
  IC50c4p_corr[index,] <- c(colnames(IC50c4p[col]),"Kendall",kendall_corr$statistic,
                        NA,kendall_corr$p.value,NA, NA,kendall_corr$estimate)
  
  index = index + 1
}

#make it print out the ones with p<0.05
for (i in 1:nrow(IC50c4p_corr)) {
  if (IC50c4p_corr[i,5] < 0.05) {
    print(IC50c4p_corr[i,c(1,2,3,5)])
  }
}

#Correlation with IC50tet ----
IC50tet <- merge(IC50tet,IC50, by="Snake")
IC50tet <- IC50tet[!is.na(IC50tet$IC50),]

#linear regressions
for (col in 4:13) {
  plot(IC50tet$IC50,IC50tet[,col],main = colnames(IC50tet)[col], xlab = "TTX Resistance (MAMU)",
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

#Correlation with rheobase ----