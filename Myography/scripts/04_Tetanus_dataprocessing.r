# TetanusP2 - Data processing
#   Protocol 2: One 2sec pulse
# 1st stimulus occurs at 0.5s (seconds)
# System delay is 4ms (miliseconds)
# Pulsewidth is 500.0us (microseconds)
# Stimulus rate is 1000 Hz
# Last stimulus is at 2.5s
#
# Processing:
# 1. Pulse on if stimulus > 1.0
# 2. Base force = average of prior 1000 readings
# 3. Data of interest - 30,000 (3 secs) from start of pulse
# 4. Normalized force = (measuredF - baseF)*9.80665/(Muscle Mass in grams)
#

# Get Snake Info - Species, Genotype etc.
sinf = read.csv("./data_raw/Snake_data_sheets/SnakeInfo-09.30.2020.csv")

# Get Snake Muscle Masses
smm = read.csv("./data_raw/Snake_data_sheets/SnakeSkeletalMuscleMasses-9.28.2020.csv")
# reset all missing muscle mass values to -1.0smm
smm[is.na(smm)] <- 0.999

# dname = "Tcouchii_side_hustle/Couchii_tetanus"
# # dname = "C:/Bobby/Data-CSV/p3Tetanus"
# setwd(dname)
# 
# x = strsplit(dname, "/")

# Set up file for output summary metrics (as we did with C4P)
ofsum = "OutFiles/Tetanus/test/Couchii_Tetanus_Metrics.csv"
hdrs <-
  c(
    "Species",
    "Snake",
    "Muscle",
    "BaseF(N/g)",
    "ContrAmpl(N/g)",
    "To10pct(ms)",
    "To50pct(ms)",
    "10to50pct(ms)",
    "FMaxRateOfChg(N/g/s)",
    "FMinRateOfChg(N/g/s)",
    "ToFChgMax(ms)",
    "ToFChgMin(ms)",
    "DiffFChgMaxToMin(ms)",
    "MAMU",
    "MusMassg"
  )
write.table(
  t(hdrs),
  file = ofsum,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

#set up output file for force data
ofF <- "OutFiles/Tetanus/test/Couchii_Tetanus_Force.csv"
ofFhdr <-
  c(
    "Species",
    "Snake",
    "Muscle",
    "MAMU",
    "MusMassg",
    format((1:30000) / 10000, scientific = F)
  )
write.table(
  t(ofFhdr),
  file = ofF,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

#Set up file for first derivative of force results
ofF1d <- "OutFiles/Tetanus/test/Couchii_Tetanus_Force_1d.csv"
ofF1dhdr <-
  c(
    "Species",
    "Snake",
    "Muscle",
    "MAMU",
    "MusMassg",
    (1:599) / 200
  )
write.table(
  t(ofF1dhdr),
  file = ofF1d,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

# Make a list of all the raw data files
files = list.files(path = "data_raw/Tetanus/", pattern = "csv")

#Set up a file to track which files we've gone through
fTet <- "OutFiles/Tetanus/test/P2-TetanusFiles.csv"
write.table(
  t(c("File", "Snake", "Muscle", "Mass(g)")),
  file = fTet,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

#Now we go through each raw data file, calculate all our metrics, and append to
# the storage files we made earlier
for (file in files) {
  y <- strsplit(file, "-")
  #Extract info from the name of the file
  
  snake <- trimws(y[[1]][1])
  # get Snake Info
  sSpecies <- as.character(sinf$Species[which(sinf$Snake == snake)])
  sMAMU <- as.character(sinf$MAMU[which(sinf$Snake == snake)])
  
  # get Muscle Mass
  muscle <- trimws(y[[1]][2])
  sMusMassg <- smm[which(smm$SnakeID == snake), muscle] / 1000
  
  write.table(
    t(c(file, snake, muscle, sMusMassg)),
    file = fTet,
    append = TRUE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE
  )
  
  #ofdata <- paste("./output/", snake, "-", muscle, ".csv", sep = "")
  #Actually read in the data file
  raw = read.csv(paste("data_raw/Tetanus/",file, sep =""))
  
  # Select all rows with stimulus > 1 (and filter out rows with insufficient stimulus)
  stmRows <- as.integer(rownames(raw[raw$Stimulus > 1.0, ]))
  if (length(stmRows) < 99) {
    cat("File: ", file, "; Stim Rows: ", length(stmRows), "\n")
    next
  }
  # extract force-baseline for 1000 pts before start of stimulus;
  sr <- stmRows[1]
  meanF = mean(raw$Force..g.[sr - 1000:sr])
  rspF <-
    (raw$Force..g.[stmRows[1]:(stmRows[1] + 29999)] - meanF) * (9.80665 / sMusMassg)
  oFline <- c(
    sSpecies,
    snake,
    muscle,
    sMAMU,
    sMusMassg,
    rspF
  )
  write.table(
    t(oFline),
    file = ofF,
    append = TRUE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE
  )
  #  Low pass 200Hz filter -> pick every 30th entry & convert to per sec
  # This will give us our first derivative file
  rspF1d <- diff(rspF[1:600 * 50]) * 200
  oF1dline <- c(
    sSpecies,
    snake,
    muscle,
    sMAMU,
    sMusMassg,
    rspF1d
  )
  write.table(
    t(oF1dline),
    file = ofF1d,
    append = TRUE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE
  )
  #Calculating all the summary stats (again, just like the C4P script)
  tmf <- max(round(rspF, 0)) - 10
  tmmin = min(which(rspF >= tmf))
  tmmax <- max(which(rspF >= tmf))
  maxF <- mean(rspF[tmmin:tmmax])
  t100p <- round((tmmin + tmmax) / 2)
  # convert all times to ms from tenth of a ms
  t10p <- min(which(rspF > 0.1 * maxF))
  t50p <- t100p + min(which(rspF[t100p:length(rspF)] <= 0.5 * maxF))
  maxF1d <- max(rspF1d)
  minF1d <- min(rspF1d)
  t1dmax <- min(which(rspF1d >= maxF1d)) * 50
  t1dmin <- min(which(rspF1d <= minF1d)) * 50
  sumLine = c(
    sSpecies,
    snake,
    muscle,
    meanF,
    maxF,
    t10p / 10,
    t50p / 10,
    (t50p - t10p) / 10,
    maxF1d,
    minF1d,
    t1dmax / 10,
    t1dmin / 10,
    (t1dmin - t1dmax) / 10,
    sMAMU,
    sMusMassg
  )
  write.table(
    t(sumLine),
    file = ofsum,
    append = TRUE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE
  )
}

#not really anything to sort here but I'm still making a flipped version to maybe 
# use for analyses (a bit slow because of the size of the file)
df = read.csv("OutFiles/Tetanus/test/Couchii_Tetanus_Force.csv")
#write.csv(t(df[order(df$Snake, df$Muscle),]),
#          "OutFiles/Tetanus/test/Couchii_Tetanus_Force_Sorted.csv")

#filter out outliers/invalid data based on waveform ----
# library(ggplot2)
# library(tidyverse)
# 
# #some muscles likely tore themselves or had other issues based on waveform, 
# # and there's few enough runs I can filter those out manually instead of writing
# # some sort of algorithm to recognize bad waveforms
# 
# #need to plot only 5 snakes at a time for colors to be distinguishable
# #(i promise i'm aware of how incredibly cursed this method is, please forgive me)
# plot_IDs_1 <- ggplot(subset(tet_force_long, 
#                             Snake %in% c("CRF2630","CRF2631","CRF2633", "CRF2669","CRF2670")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_1
# #lots of bad ones, 2669 and 2631 both look ripped and 2670 is a bit small
# 
# plot_IDs_2 <- ggplot(subset(tet_force_long, 
#                             Snake %in% c("CRF2671", "CRF2672", "CRF2673", "CRF2674", "CRF2675")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_2
# #all look fine, though the normalization is horrid
# #2675 is the tall outlier, contrary to bobby's suggestion it's actually on the smaller side (11 mg)
# 
# plot_IDs_3 <- ggplot(subset(tet_force_long, 
#                             Snake %in% c("CRF2676", "CRF2677", "CRF2678", "CRF2679", "CRF2680")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_3
# #CRF2680 ripped
# 
# plot_IDs_4 <- ggplot(subset(tet_force_long, 
#                             Snake %in% c("CRF2681", "CRF3051", "CRF3052", "CRF3055", "CRF3056")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_4
# #all fine (except maybe 2681)
# 
# plot_IDs_5 <- ggplot(subset(tet_force_long, 
#                             Snake %in% c("CRF3065", "CRF3066", "CRF3070", "CRF3074")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_5
# #3074 totally wrong, 3065 and 3066 also look bad
# 
# plot_IDs_6 <- ggplot(subset(tet_force_long, 
#                             Snake %in% c("CRF3058", "CRF3059", "CRF3060", "CRF3061", "CRF3064")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_6
# #3060 has weird decay
# tet_force_long_sub <- tet_force_long %>% 
#   filter(!Snake %in% c("CRF3060","CRF3074","CRF3065","CRF3066","CRF2680","CRF2669","CRF2631","CRF2670")) 
# #%>% filter(Snake != "CRF2675")
# 
# #save this file
# #write.csv(tet_force_long_sub, "data_processed/tetanus_no_outliers.csv")