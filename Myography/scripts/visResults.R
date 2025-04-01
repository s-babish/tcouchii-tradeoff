#basic visualizing results before actual statistical analysis (mostly just scatterplots)

setwd("C:/Users/sdbab/OneDrive - University of Nevada, Reno/UNR/tcouchii/Myography")

#C4P Results ----
c4p <- read.csv("Tcouchii_side_hustle/Couchii_C4P/output/Couchii_C4P.csv")
head(c4p)

for (col in 7:17) {
  plot(c4p$MAMU,c4p[,col],main = colnames(c4p)[col], xlab = "TTX Resistance (MAMU)",
       ylab = colnames(c4p)[col])
  model <- lm(c4p[,col] ~ c4p$MAMU)
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
#highest r2 is 0.03 

#Tetanus results ----
tetanus <- read.csv("Tcouchii_side_hustle/Couchii_tetanus/output/Couchii_tetanus.csv")
head(tetanus)

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
#highest r2 is 0.05

#Rheobase results (kind of) -----

pulse_l <- c(50,100,200,500,1000,10000,50000)

for (i in 1:length(pulse_l)) {
  infile <- toString(paste(pulse_l[i],"_TPA-sigmoidal-rpt.csv", sep=""))
  df <- read.csv(infile)
  par(mfrow=c(1,3))
  for (col in 3:5) {
    plot(df$MAMU,df[,col],main = paste(colnames(df)[col]," for pulse length ",
                                       pulse_l[i]), xlab = "TTX Resistance (MAMU)",
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
    legend(x = "topright", bty = "n",
           legend = bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))) 
    }
}

#highest r2 is 0.14