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

#original -----
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
                      control = nls.lm.control(maxiter = 10000, maxfev = 10000), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
    print(nls.out)
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

#transform x0 via p0 ----
for (i in 1:length(pulse_l)) {
  
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile, row.names = 2)
  df <- df[1:nrow(df),2:ncol(df)]
  
  sigmoidal <- function(x,A1,A2,x0,dx){
    (A1-A2)/(1 + exp((log(x)-log(x0))/dx)) + A2
  }
  sigmoidal <- Vectorize(sigmoidal)
  
  varsigmoidal <- function(x,p0,dx){
    x0 <- exp(p0)
    out <- (1)/(1 + exp((log(x)-log(x0))/dx))
    return(out)
  }
  varsigmoidal <- Vectorize(varsigmoidal)
  
  residFun <- function(parS,observed,indices){
    sigmoidal(as.numeric(rownames(df))[indices],parS$A1,parS$A2,parS$x0,parS$dx) - observed
  }
  
  varresidFun <- function(parS,observed,indices){
    varsigmoidal(as.numeric(rownames(df))[indices],parS$p0,parS$dx) - observed
  }
  
  parStart <- list(A2 = 1, A1 = 0, p0 = 5 , dx = 4) #as dx-->0, steepness --> infinity; was x0 = 60
  fitParams <- function(x){
    nls.out <- nls.lm(par = parStart, fn = varresidFun, 
                      control = nls.lm.control(maxiter = 10000, maxfev = 10000), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
    print(nls.out)
    unlist(nls.out$par[3:4])
  }
  
  result <- Vectorize(fitParams)(1:ncol(df))
  colnames(result) <- colnames(df)
  print(result)
}

#modify residuals/likelihood function (also I realized nls.lm allows bounds on params) -----
for (i in 1:length(pulse_l)) {
  
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile, row.names = 2)
  df <- df[1:nrow(df),2:ncol(df)]
  
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
    if(parS$x0 > 0) {
      return(varsigmoidal(as.numeric(rownames(df))[indices],parS$x0,parS$dx) - observed)
    } else {
      return(100000000)
    }
  }
  
  parStart <- list(A2 = 1, A1 = 0, x0 = 60 , dx = 4) #as dx-->0, steepness --> infinity
  fitParams <- function(x){
    nls.out <- nls.lm(par = parStart, fn = varresidFun, 
                      control = list(maxiter = 10000, maxfev = 10000),
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]),
                      lower = c(0,0,0,0), upper = c(10,10,600,10000000000))
    print(nls.out)
    plot(1:(nls.out$niter+1), log(nls.out$rsstrace), type = "b")
    unlist(nls.out$par[3:4])
  }
  
  result <- Vectorize(fitParams)(1:ncol(df))
  colnames(result) <- colnames(df)
  print(result)
}

#optim (last ditch resort) -----
for (i in 1:length(pulse_l)) {
  
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile, row.names = 2)
  df <- df[1:nrow(df),2:ncol(df)]
  rownames(df) <- as.numeric(rownames(df))
  rownames(df)[rownames(df) == 0] <- 1e-6  # Replace 0 with a small number
  
  
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
  
  varresidFun <- function(par, observed, indices){
    x0 <- par[3]
    dx <- par[4]
    x_vals <- as.numeric(rownames(df))[indices]
    preds <- varsigmoidal(x_vals, x0, dx)
    sum((preds - observed)^2)  # Return sum of squared residuals
  }
  
  
  parStart <- c(A2 = 1, A1 = 0, x0 = 60 , dx = 4) #as dx-->0, steepness --> infinity
  fitParams <- function(x){
    observed <- df[,x]
    indices <- which(!is.na(observed))
    observed <- observed[indices]
    
    optim.out <- optim(par = parStart, 
                       fn = varresidFun, 
                       observed = observed, 
                       indices = indices,
                      control = list(maxit = 10000, trace = 3),
                      method = "L-BFGS-B",
                      lower = c(0,0,0,0), 
                      upper = c(10,10,600,Inf))
    print(optim.out)
    unlist(optim.out$par[3:4])
  }
  result <- sapply(1:ncol(df), fitParams)
  colnames(result) <- colnames(df)
  print(result)

}

#I think either something is going very wrong in the process or these data just 
# do not fit this equation the way Bobby wants them too
#will do some research to see if there are other commonly used equations out there

#fixed the parameters -----
for (i in 1:length(pulse_l)) {
  
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile, row.names = 2)
  df <- df[1:nrow(df),2:ncol(df)]
  
  #realized in the process of doing this that the reason this hasn't been working 
  # is that the parameters are backwards! I hate it here!
  
  #A1 is min, A2 is max
  sigmoidal <- function(x,A1,A2,x0,dx){
    (A1-A2)/(1 + exp((log(x)-log(x0))/dx)) + A2
  }
  sigmoidal <- Vectorize(sigmoidal)
  
  varsigmoidal <- function(x,x0,dx){
    (-1)/(1 + exp((log(x)-log(x0))/dx)) + 1
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
                      control = nls.lm.control(maxiter = 1024, maxfev = 1024), 
                      observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]),
                      lower = c(0,0,0,NA), upper = c(10,10,600,NA))
    print(nls.out)
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
  maxslope <- apply(result,2,function(p){max.abs(sigmoidalDeriv(as.numeric(rownames(df)),0,1,p[1],p[2]), na.rm = T)})
  
  solveSigmoidal <- function(y,A1,A2,x0,dx){
    exp(log(-(A2 - A1)/(y - A2) - 1) * dx + log(x0))
  }
  solveSigmoidal <- Vectorize(solveSigmoidal)
  range.10.90 <- abs(solveSigmoidal(0.9,0,1,result[1,],result[2,]) - solveSigmoidal(0.1,0,1,result[1,],result[2,]))
  
  result <- rbind(result,maxslope,range.10.90)
  
  filename <- toString(paste("OutFiles/Rheobase/Sigmoidal/test/",pulse_l[i],
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
  
  filename2 <- toString(paste("OutFiles/Rheobase/Sigmoidal/test/",pulse_l[i],
                              "_Sigmoidal-rpt.csv", sep=""))
  write.csv(rdf, filename2,row.names = F)
}

#I realized this doesn't make the normal "TPA-sigmoidal-rpt.csv" file, I need
# to do that
#actually i was just using that for snake ID and MAMU so I can join that another way

rheobase_sigmoidal <- matrix(ncol=ncol(rdf)+1)
colnames(rheobase_sigmoidal)<-c(colnames(rdf),"pulse_length")

for (i in 1:length(pulse_l)) {
  infile <- toString(paste("OutFiles/Rheobase/Sigmoidal/test/",
                           pulse_l[i],"_Sigmoidal-rpt.csv", sep=""))
  df <- read.csv(infile)
  df_pulseID <- df %>% 
    mutate(
      pulse_length = pulse_l[i] #add column to differentiate entries by pulse length
    )
  
  #join to master storage matrix
  rheobase_sigmoidal <- rbind(rheobase_sigmoidal,df_pulseID)
}

write.csv(rheobase_sigmoidal, "OutFiles/Rheobase/Sigmoidal/test/Sigmoidal-rpt.csv")
