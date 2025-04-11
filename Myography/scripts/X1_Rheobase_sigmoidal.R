#Libraries -----
library(tidyverse)
library(minpack.lm)

#Load (some) data ----
# Get Snake Info
sinf = read.csv("./data_raw/Snake_data_sheets/SnakeInfo-09.30.2020.csv")

# Get Snake Muscle Masses
smm = read.csv("./data_raw/Snake_data_sheets/SnakeSkeletalMuscleMasses-9.28.2020.csv")
smm[is.na(smm)] <- 0.999 # reset all missing muscle mass values to -1.0smm

#Compile and format data we actually want
smm_gathered <- smm %>% gather(Muscle, MusMassg, M1:M13) %>% select(-Date) %>% rename(Snake = SnakeID)


pulse_l <- c(50,100,200,500,1000,10000,50000)
for (i in 1:length(pulse_l)) {
 
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile, row.names = 2)
  df <- df[1:nrow(df),2:ncol(df)]
  
  #bobby had this chunk commented out already so I'm not worried about it
  #df <- df[,which(colSums(!is.na(df)) >= 4)] #it appears to be asking how many NAs are in each column, expecting great than or equal to 4?
  #the base function which() will return a value for the position of something (e.g. Which(letters == g) returns 7)

  sigmoidal <- function(x,A1,A2,x0,dx){
    (A1-A2)/(1 + exp((log(x)-log(x0))/dx)) + A2
  }
  sigmoidal <- Vectorize(sigmoidal)
  
  varsigmoidal <- function(x,x0,dx){
    (1)/(1 + exp((log(x)-log(x0))/dx))
  }
  varsigmoidal <- Vectorize(varsigmoidal)
  
  residFun <- function(parS,observed,indices){
    sigmoidal(as.numeric(rownames(df))[indices],parS$A1,parS$A2,parS$x0,parS$dx) - observed
  }
  
  varresidFun <- function(parS,observed,indices){
    varsigmoidal(as.numeric(rownames(df))[indices],parS$x0,parS$dx) - observed
  }
  
  parStart <- list(A2 = 1, A1 = 0, x0 = 60 , dx = 4) #as dx-->0, steepness --> infinity
  fitParams <- function(x){
    nls.out <- nls.lm(par = parStart, fn = varresidFun, 
                      control = nls.lm.control(maxiter = 100), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
    unlist(nls.out$par[3:4])
  }
  
  result <- Vectorize(fitParams)(1:ncol(df))
  colnames(result) <- colnames(df)
  print(result)
  sigmoidalDeriv <- function(x,A1,A2,x0,dx){
    -((A1 - A2) * exp((log(x) + log(x0))/dx))/(dx * x * (exp(log(x)/dx) + exp(log(x0)/dx))^2)
  }
  sigmoidalDeriv <- Vectorize(sigmoidalDeriv)
  
  max.abs <- function(x,...){
    sign(x[which.max(abs(x))])*max(abs(x),...)
  }
  maxslope <- apply(result,2,function(p){max.abs(sigmoidalDeriv(as.numeric(rownames(df)),1,0,p[1],p[2]), na.rm = T)})
  
  solveSigmoidal <- function(y,A1,A2,x0,dx){
    exp(log(-(A2 - A1)/(y - A2) - 1) * dx + log(x0))
  }
  solveSigmoidal <- Vectorize(solveSigmoidal)
  range.10.90 <- abs(solveSigmoidal(0.9,1,0,result[1,],result[2,]) - solveSigmoidal(0.1,1,0,result[1,],result[2,]))
  
  result <- rbind(result,maxslope,range.10.90)
  
  filename <- toString(paste(pulse_l[i],"_TPA-sigmoidal.csv", sep=""))
  #write.csv(t(result), filename)
  
  rdf <- read.csv(filename)
  x <- strsplit(as.character(rdf$X),"_")
  rdf$Snake <- unlist(lapply(x,'[[',1))
  rdf$Muscle <- unlist(lapply(x,'[[',2))
  #rdf$Temperature <- unlist(lapply(x,'[[',3)) no temp data ----
  rdf$Species <- sinf$Species[match(rdf$Snake, sinf$Snake)]
  rdf$Genotype <- sinf$Genotype[match(rdf$Snake, sinf$Snake)]
  rdf$MAMU <- sinf$MAMU[match(rdf$Snake, sinf$Snake)]
  rdf$County <- sinf$COUNTY[match(rdf$Snake, sinf$Snake)]
  rdf$Long <- sinf$Longitude[match(rdf$Snake, sinf$Snake)]
  rdf$Lat <- sinf$Latitude[match(rdf$Snake, sinf$Snake)]
  rdf$SVLmm <- as.character(sinf$SVLmm[match(rdf$Snake, sinf$Snake)])
  rdf <- rdf %>% left_join(smm_gathered)
  rdf$BodMassg <- as.character(sinf$BodMassg[match(rdf$Snake, sinf$Snake)])
  rdf$Sex <- as.character(sinf$Sex[match(rdf$Snake, sinf$Snake)])
  rdf$DtExp <- as.character(sinf$Date_Experimented[match(rdf$Snake, sinf$Snake)])
  rdf$DtColl <- as.character(sinf$Date_Collected[match(rdf$Snake, sinf$Snake)])
  rdf$DInBet <- as.character(sinf$Days_in_Between[match(rdf$Snake, sinf$Snake)])
  
  filename2 <- toString(paste(pulse_l[i],"_TPA-sigmoidal-rpt.csv", sep=""))
  #write.csv(rdf, filename2)
}
