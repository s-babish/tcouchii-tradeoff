#Rheobase correlation results (kind of) -----
#doesn't work rn
# pulse_l <- c(50,100,200,500,1000,10000,50000)
# 
# for (i in 1:length(pulse_l)) {
#   infile <- toString(paste(pulse_l[i],"_TPA-sigmoidal-rpt.csv", sep=""))
#   df <- read.csv(infile)
#   par(mfrow=c(1,3))
#   for (col in 3:5) {
#     print(pulse_l, colnames(df[col]))
#     
#     pearson_corr <- cor.test(df$MAMU,df[,col], method = "pearson")
#     print(pearson_corr)
#     
#     kendall_corr <- cor.test(df$MAMU,df[,col], method = "kendall")
#     print(kendall_corr)
#   }
# }
# 

#Read in and compile rheobase sigmoidal results (calculated by script TPA_sigmoidal_2024.R)
pulse_l <- c(50,100,200,500,1000,10000,50000)
i=1

keydf <- read.csv("TPA-sigmoidal-rpt.csv")
keydf <- keydf[,c(7,11)]

#create a storage matrix to hold all results in
rheobase_sigmoidal <- matrix(ncol=6)
colnames(rheobase_sigmoidal)<-c("X","x0","dx","maxslope","range.10.90","pulse_length")

for (i in 1:length(pulse_l)) {
  infile <- toString(paste(pulse_l[i],"_TPA-sigmoidal.csv", sep=""))
  df <- read.csv(infile)
  df_pulseID <- df %>% 
    mutate(
      pulse_length = pulse_l[i] #add column to differentiate entries by pulse length
    )
  
  #join to master storage matrix
  rheobase_sigmoidal <- rbind(rheobase_sigmoidal,df_pulseID)
}

rheobase_sigmoidal <- rheobase_sigmoidal %>% 
  drop_na() %>% 
  mutate(
        Snake = sub("_[^_]+$", "", X),
        MAMU = 0
      ) %>%
  rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
