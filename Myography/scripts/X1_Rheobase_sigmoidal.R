#Libraries -----
library(tidyverse)
library(minpack.lm)

#Load (some) data ----
#give folder name (change according to dataset you're using)
foldername <- "P"

# Get Snake Info
sinf = read.csv("./data_raw/Snake_data_sheets/SnakeInfo-09.30.2020.csv")

# Get Snake Muscle Masses
smm = read.csv("./data_raw/Snake_data_sheets/SnakeSkeletalMuscleMasses-9.28.2020.csv")
smm[is.na(smm)] <- 0.999 # reset all missing muscle mass values to -1.0smm

#Compile and format data we actually want
smm_gathered <- smm %>% gather(Muscle, MusMassg, M1:M13) %>% 
  #select(-Date) %>% 
  rename(Snake = SnakeID)


pulse_l <- c(50,100,200,500,1000,10000,50000)


#fixed the parameters -----
for (i in 1:length(pulse_l)) {
  
  infile <- toString(paste("OutFiles/Rheobase/",foldername,"/Force_scaled/Rheo_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile, row.names = 2)
  
  #save max values and then drop them from the fitting
  maxes <- df[17,-1]
  df <- df[1:(nrow(df)-1),2:ncol(df)]
  
  
  #realized in the process of doing this that the reason this hasn't been working 
  # is that the parameters are backwards! I hate it here!
  
  #and the dx is fucky and written such that higher dx is actually a flatter curve, 
  #not a steeper one! I'm changing the equation I use in the variant to the one
  #most people actually use but leaving the old one for posterity
  
  #(note: the old one was a log-logistic regression, whereas the new one is a 
  #generalized logistic/hill form)
  
  #A1 is min, A2 is max
  # sigmoidal <- function(x,A1,A2,x0,dx){
  #   (A1-A2)/(1 + exp((log(x)-log(x0))/dx)) + A2
  # }
  # sigmoidal <- Vectorize(sigmoidal)
  # 
  # varsigmoidal <- function(x,x0,dx){
  #   (-1)/(1 + exp((log(x)-log(x0))/dx)) + 1
  # }
  
  varsigmoidal <- function(x,x0,dx) {
    1 + (-1) / (1 + (x / x0)^dx)
  }
  varsigmoidal <- Vectorize(varsigmoidal)
  
  varresidFun <- function(parS,observed,indices){
    varsigmoidal(as.numeric(rownames(df))[indices],parS$x0,parS$dx) - observed
  }
  
  parStart <- list(A2 = 1, A1 = 0, x0 = 60 , dx = 4) 
  fiRheorams <- function(x){
    nls.out <- nls.lm(par = parStart, fn = varresidFun, 
                      control = nls.lm.control(maxiter = 1024, maxfev = 1024), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]),
                      lower = c(0,0,0,NA), upper = c(10,10,600,NA))
    print(nls.out)
    unlist(nls.out$par[3:4])
  }
  
  result <- Vectorize(fiRheorams)(1:ncol(df))
  colnames(result) <- colnames(df)
  print(result)
  
  #also changing the old version of this
  # sigmoidalDeriv <- function(x,A1,A2,x0,dx){
  #   -((A1 - A2) * exp((log(x) + log(x0))/dx))/(dx * x * (exp(log(x)/dx) + exp(log(x0)/dx))^2)
  # }
  
  sigmoidalDeriv <- function(x, A1, A2, x0, dx) {
    numerator <- (A1 - A2) * dx * x^(dx - 1)
    denominator <- x0^dx * (1 + (x / x0)^dx)^2
    -numerator / denominator
  }
  
  sigmoidalDeriv <- Vectorize(sigmoidalDeriv)
  
  max.abs <- function(x,...){
    sign(x[which.max(abs(x))])*max(abs(x),...)
  }
  maxslope <- apply(result,2,function(p){max.abs(sigmoidalDeriv(as.numeric(rownames(df)),0,1,p[1],p[2]), na.rm = T)})
  
  #and this one
  # solveSigmoidal <- function(y,A1,A2,x0,dx){
  #   exp(log(-(A2 - A1)/(y - A2) - 1) * dx + log(x0))
  # }
  
  solveSigmoidal <- function(y,A1,A2,x0,dx){
    x0 * ((A1 - y) / (y - A2))^(1 / dx)
  }
  solveSigmoidal <- Vectorize(solveSigmoidal)
  ec10 <- solveSigmoidal(0.1,0,1,result[1,],result[2,]) #adding at bobby's suggestion
  ec90 <- solveSigmoidal(0.9,0,1,result[1,],result[2,]) #ibid
  range.10.90 <- abs(ec90 - ec10)
  
  result <- rbind(result,maxslope,ec10,ec90,range.10.90,maxes)
  rownames(result) <- c("x0","dx","maxslope","EC10","EC90","range1090","maxF")
  
  filename <- toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/",pulse_l[i],
                             "_Sigmoidal.csv", sep=""))
  write.csv(t(result), filename)
  
  rdf <- read.csv(filename)
  x <- strsplit(as.character(rdf$X),"_")
  rdf$Snake <- unlist(lapply(x,'[[',1))
  rdf$Muscle <- unlist(lapply(x,'[[',2))
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
  
  filename2 <- toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/",pulse_l[i],
                              "_Sigmoidal-rpt.csv", sep=""))
  write.csv(rdf, filename2,row.names = F)
}

#make general sigmoidal report file ----

rheobase_sigmoidal <- matrix(ncol=ncol(rdf)+1)
colnames(rheobase_sigmoidal)<-c(colnames(rdf),"pulse_length")

for (i in 1:length(pulse_l)) {
  infile <- toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/",
                           pulse_l[i],"_Sigmoidal-rpt.csv", sep=""))
  df <- read.csv(infile)
  df_pulseID <- df %>% 
    mutate(
      pulse_length = pulse_l[i] #add column to differentiate entries by pulse length
    )
  
  #join to master storage matrix
  rheobase_sigmoidal <- rbind(rheobase_sigmoidal,df_pulseID)
}

outfile <- toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv", sep=""))
write.csv(rheobase_sigmoidal, outfile)

#log the current values (x) ----
# for (i in 1:length(pulse_l)) {
#   
#   infile <- toString(paste("OutFiles/Rheobase/WT/Rheo_",pulse_l[i],".csv", sep=""))
#   df <- read.csv(infile, row.names = 2)
#   
#   #save max values and then drop them from the fitting
#   maxes <- df[17,-1]
#   df <- df[1:(nrow(df)-1),2:ncol(df)]
#   
#   rownames(df) <- c(1,log(as.numeric(rownames(df)[-1])))
#   
#   #bobby suggested logging the current which might have been what the old eq was 
#   # but in a weirder way; anyway i'm trying it
#   
#   #can't have a 0 in a log
#   varsigmoidal <- function(x,x0,dx) {
#     1 + (-1) / (1 + (x / x0)^dx)
#   }
#   varsigmoidal <- Vectorize(varsigmoidal)
#   
#   varresidFun <- function(parS,observed,indices){
#     varsigmoidal(as.numeric(rownames(df))[indices],parS$x0,parS$dx) - observed
#   }
#   
#   parStart <- list(A2 = 1, A1 = 0, x0 = 4 , dx = 4) 
#   fiRheorams <- function(x){
#     nls.out <- nls.lm(par = parStart, fn = varresidFun, 
#                       control = nls.lm.control(maxiter = 1024, maxfev = 1024), 
#                       observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]),
#                       lower = c(0,0,0,NA), upper = c(10,10,600,NA))
#     print(nls.out)
#     unlist(nls.out$par[3:4])
#   }
#   
#   result <- Vectorize(fiRheorams)(1:ncol(df))
#   colnames(result) <- colnames(df)
#   print(result)
#   
#   #this doesn't make sense rn but i'll fix it later if it looks like this is gonna work
#   
#   sigmoidalDeriv <- function(x, A1, A2, x0, dx) {
#     numerator <- (A1 - A2) * dx * x^(dx - 1)
#     denominator <- x0^dx * (1 + (x / x0)^dx)^2
#     -numerator / denominator
#   }
#   
#   sigmoidalDeriv <- Vectorize(sigmoidalDeriv)
#   
#   max.abs <- function(x,...){
#     sign(x[which.max(abs(x))])*max(abs(x),...)
#   }
#   maxslope <- apply(result,2,function(p){max.abs(sigmoidalDeriv(as.numeric(rownames(df)),0,1,p[1],p[2]), na.rm = T)})
#   
#   
#   solveSigmoidal <- function(y,A1,A2,x0,dx){
#     exp(x0 * ((A1 - y) / (y - A2))^(1 / dx))
#   }
#   solveSigmoidal <- Vectorize(solveSigmoidal)
#   range.10.90 <- abs(solveSigmoidal(0.9,0,1,result[1,],result[2,]) - solveSigmoidal(0.1,0,1,result[1,],result[2,]))
#   
#   result <- rbind(result,maxslope,range.10.90,maxes)
#   rownames(result) <- c("x0","dx","maxslope","range1090","maxF")
#   
#   filename <- toString(paste("OutFiles/Rheobase/Sigmoidal/logtest/",pulse_l[i],
#                              "_Sigmoidal.csv", sep=""))
#   write.csv(t(result), filename)
#   
#   rdf <- read.csv(filename)
#   x <- strsplit(as.character(rdf$X),"_")
#   rdf$Snake <- unlist(lapply(x,'[[',1))
#   rdf$Muscle <- unlist(lapply(x,'[[',2))
#   rdf$Species <- sinf$Species[match(rdf$Snake, sinf$Snake)]
#   rdf$Genotype <- sinf$Genotype[match(rdf$Snake, sinf$Snake)]
#   rdf$MAMU <- sinf$MAMU[match(rdf$Snake, sinf$Snake)]
#   rdf$County <- sinf$COUNTY[match(rdf$Snake, sinf$Snake)]
#   rdf$Long <- sinf$Longitude[match(rdf$Snake, sinf$Snake)]
#   rdf$Lat <- sinf$Latitude[match(rdf$Snake, sinf$Snake)]
#   rdf$SVLmm <- as.character(sinf$SVLmm[match(rdf$Snake, sinf$Snake)])
#   rdf <- rdf %>% left_join(smm_gathered)
#   rdf$BodMassg <- as.character(sinf$BodMassg[match(rdf$Snake, sinf$Snake)])
#   rdf$Sex <- as.character(sinf$Sex[match(rdf$Snake, sinf$Snake)])
#   rdf$DtExp <- as.character(sinf$Date_Experimented[match(rdf$Snake, sinf$Snake)])
#   rdf$DtColl <- as.character(sinf$Date_Collected[match(rdf$Snake, sinf$Snake)])
#   rdf$DInBet <- as.character(sinf$Days_in_Between[match(rdf$Snake, sinf$Snake)])
#   
#   filename2 <- toString(paste("OutFiles/Rheobase/Sigmoidal/logtest/",pulse_l[i],
#                               "_Sigmoidal-rpt.csv", sep=""))
#   write.csv(rdf, filename2,row.names = F)
# }


# #original -----
# for (i in 1:length(pulse_l)) {
#  
#   infile <- toString(paste("OutFiles/Rheobase/Force_scaled/test/Rheo_",pulse_l[i],".csv", sep=""))
#   df <- read.csv(infile, row.names = 2)
#   df <- df[1:nrow(df),2:ncol(df)]
#   
#   #bobby had this chunk commented out already so I'm not worried about it
#   #df <- df[,which(colSums(!is.na(df)) >= 4)] #it appears to be asking how many NAs are in each column, expecting great than or equal to 4?
#   #the base function which() will return a value for the position of something (e.g. Which(letters == g) returns 7)
# 
#   sigmoidal <- function(x,A1,A2,x0,dx){
#     (A1-A2)/(1 + exp((log(x)-log(x0))/dx)) + A2
#   }
#   sigmoidal <- Vectorize(sigmoidal)
#   
#   varsigmoidal <- function(x,x0,dx){
#     (1)/(1 + exp((log(x)-log(x0))/dx))
#   }
#   varsigmoidal <- Vectorize(varsigmoidal)
# 
#   residFun <- function(parS,observed,indices){
#     sigmoidal(as.numeric(rownames(df))[indices],parS$A1,parS$A2,parS$x0,parS$dx) - observed
#   }
#   
#   varresidFun <- function(parS,observed,indices){
#     varsigmoidal(as.numeric(rownames(df))[indices],parS$x0,parS$dx) - observed
#   }
#   
#   parStart <- list(A2 = 1, A1 = 0, x0 = 60 , dx = 4) #as dx-->0, steepness --> infinity
#   fiRheorams <- function(x){
#     nls.out <- nls.lm(par = parStart, fn = varresidFun, 
#                       control = nls.lm.control(maxiter = 10000, maxfev = 10000), 
#                       observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
#     print(nls.out)
#     unlist(nls.out$par[3:4])
#   }
#   
#   result <- Vectorize(fiRheorams)(1:ncol(df))
#   colnames(result) <- colnames(df)
#   print(result)
#   sigmoidalDeriv <- function(x,A1,A2,x0,dx){
#     -((A1 - A2) * exp((log(x) + log(x0))/dx))/(dx * x * (exp(log(x)/dx) + exp(log(x0)/dx))^2)
#   }
#   sigmoidalDeriv <- Vectorize(sigmoidalDeriv)
# 
#   max.abs <- function(x,...){
#     sign(x[which.max(abs(x))])*max(abs(x),...)
#   }
#   maxslope <- apply(result,2,function(p){max.abs(sigmoidalDeriv(as.numeric(rownames(df)),1,0,p[1],p[2]), na.rm = T)})
# 
#   solveSigmoidal <- function(y,A1,A2,x0,dx){
#     exp(log(-(A2 - A1)/(y - A2) - 1) * dx + log(x0))
#   }
#   solveSigmoidal <- Vectorize(solveSigmoidal)
#   range.10.90 <- abs(solveSigmoidal(0.9,1,0,result[1,],result[2,]) - solveSigmoidal(0.1,1,0,result[1,],result[2,]))
# 
#   result <- rbind(result,maxslope,range.10.90)
# 
#   filename <- toString(paste(pulse_l[i],"_Rheo-sigmoidal.csv", sep=""))
#   #write.csv(t(result), filename)
# 
#   rdf <- read.csv(filename)
#   x <- strsplit(as.character(rdf$X),"_")
#   rdf$Snake <- unlist(lapply(x,'[[',1))
#   rdf$Muscle <- unlist(lapply(x,'[[',2))
#   rdf$Species <- sinf$Species[match(rdf$Snake, sinf$Snake)]
#   rdf$Genotype <- sinf$Genotype[match(rdf$Snake, sinf$Snake)]
#   rdf$MAMU <- sinf$MAMU[match(rdf$Snake, sinf$Snake)]
#   rdf$County <- sinf$COUNTY[match(rdf$Snake, sinf$Snake)]
#   rdf$Long <- sinf$Longitude[match(rdf$Snake, sinf$Snake)]
#   rdf$Lat <- sinf$Latitude[match(rdf$Snake, sinf$Snake)]
#   rdf$SVLmm <- as.character(sinf$SVLmm[match(rdf$Snake, sinf$Snake)])
#   rdf <- rdf %>% left_join(smm_gathered)
#   rdf$BodMassg <- as.character(sinf$BodMassg[match(rdf$Snake, sinf$Snake)])
#   rdf$Sex <- as.character(sinf$Sex[match(rdf$Snake, sinf$Snake)])
#   rdf$DtExp <- as.character(sinf$Date_Experimented[match(rdf$Snake, sinf$Snake)])
#   rdf$DtColl <- as.character(sinf$Date_Collected[match(rdf$Snake, sinf$Snake)])
#   rdf$DInBet <- as.character(sinf$Days_in_Between[match(rdf$Snake, sinf$Snake)])
# 
#   filename2 <- toString(paste(pulse_l[i],"_Rheo-sigmoidal-rpt.csv", sep=""))
#   #write.csv(rdf, filename2)
# }
