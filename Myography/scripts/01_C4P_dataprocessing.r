# C4P - Data processing
#   Protocol 1: Control Four Pulse (C4P) Experiment
#   1st stimulus occurs at 5s (seconds)
#   System delay is 4ms (miliseconds)
#   Pulsewidth is 500.0us (microseconds)
#   Stimulus rate is 0.2 Hz
#   Last stimulus is at 24s
#
# Processing:
# 1. Pulse on if stimulus > 1.0
# 2. Base force = average of prior 1000 readings
# 3. Data of interest - 1500 (150ms) after end of pulse
# 4. Normalized force = (measuredF - baseF)*9.80665/(Muscle Mass in grams)
#

#Confirm working directory is main Myography folder (so we don't need to set it)
getwd()

# Get Snake Info - Species, Genotype etc.
sinf = read.csv("./data_raw/Snake_data_sheets/SnakeInfo-09.30.2020.csv")

# Get Snake Muscle Masses
smm = read.csv("./data_raw/Snake_data_sheets/SnakeSkeletalMuscleMasses-9.28.2020.csv")
# reset all missing muscle mass values to -1.0smm
smm[is.na(smm)] <- 0.999

# Set up file to store summary information about metrics we calculate for each
# pulse later

ofsum = "OutFiles/C4P/test/Couchii_C4P_Metrics.csv"
hdrs <- c(
  "Species",
  "Snake",
  "Muscle",
  "Rater",
  "MusMassg",
  "Pulse",
  "BaseF(N/g)",
  "ContrAmpl",
  "ToMaxF(ms)",
  "To10pct(ms)",
  "To50pct(ms)",
  "10to50pct(ms)",
  "FMaxRateOfChg(N/g/s)",
  "FMinRateOfChg(N/g/s)",
  "ToFChgMax(ms)",
  "ToFChgMin(ms)",
  "DiffFChgMaxToMin(ms)",
  "MAMU"
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
# Set up file to store force outputs
ofF <- "OutFiles/C4P/test/Couchii_C4P_Force.csv"
ofFhdr <- c(
  "Species",
  "MAMU",
  "MusMassg",
  "Pulse",
  "Snake",
  "Muscle",
  "Rater",
  format((1:1500) / 10000,scientific = F)
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
# And make a file to store information on the first derivative of the force
ofF1d <- "OutFiles/C4P/test/Couchii_C4P_Force_1d.csv"
ofF1dhdr <- c(
  "Species",
  "MAMU",
  "MusMassg",
  "Pulse",
  "Snake",
  "Muscle",
  "Rater",
  (1:29) / 200
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

# get all data files (in data_raw folder)
files = list.files(path = "data_raw/C4P/", pattern = "csv")

#Set up a file just to track which snakes we've gone through
fC4P <- "OutFiles/C4P/test/P1-C4PFiles.csv"
write.table(
  t(c("File", "Snake", "Muscle", "Rater", "Mass(g)")),
  file = fC4P,
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
  #extracting info from the name of the file
  y <- strsplit(file, "-")
  snake <- trimws(y[[1]][1])
  # get Snake Info
  sSpecies <- as.character(sinf$Species[which(sinf$Snake == snake)])
  sMAMU <- as.character(sinf$MAMU[which(sinf$Snake == snake)])
  # Get Muscle Mass
  muscle <- trimws(y[[1]][2])
  # Get Rater from Snake Info file for snake/musc combo
  rater <- as.character(sinf[which(sinf$Snake == snake), paste0("Rater..", muscle, ".")])

  MusMassg <- smm[which(smm$SnakeID == snake), muscle] / 1000
  write.table(
    t(c(file, snake, muscle, rater, MusMassg)),
    file = fC4P,
    append = TRUE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE
  )
  
  #actually read in the data file
  raw = read.csv(paste("data_raw/C4P/",file, sep =""))
  # Select all rows with stimulus > 1 (and filter out files without all 4 pulses)
  stmRows <- as.integer(rownames(raw[raw$Stimulus > 1.0,]))
  if (length(stmRows) < 4) {
    cat("File: ", file, "; Stim Rows: ", length(stmRows), "\n")
    next
  }
  categories <- cutree(hclust(dist(stmRows)), k = 4)
  stmRows <- aggregate(stmRows, list(cats = categories), "max")[, 2]
  if (length(stmRows) < 4) {
    cat("File: ", file, "; Stim Rows: ", length(stmRows), "\n")
    next
  }
  
  # extract force-baseline for 1500 pts after end of stimulus;
  i <- 0
  for (sr in stmRows) {
    i <- i + 1
    meanF = mean(raw$Force..g.[sr - 1000:sr])
    #calculate force, scaled by muscle mass
    rspF <-
      (raw$Force..g.[stmRows[1]:(stmRows[1] + 1499)] - meanF) * (0.00980665 /
                                                                   MusMassg)
    oFline <- c(
      sSpecies,
      sMAMU,
      MusMassg,
      i,
      snake,
      muscle,
      rater,
      rspF
    )
    #save calculated force (scaled and normalized to baseline) to file
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
    #  Low pass 200Hz filter -> pick every 50th entry and convert to Force / Sec
    rspF1d <- diff(rspF[1:30 * 50]) * 200
    oF1dline <- c(
      sSpecies,
      sMAMU,
      MusMassg,
      i,
      snake,
      muscle,
      rater,
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
    #Calculate all other metrics
    #max force:
    maxF <- max(rspF)
    #time to max force:
    t100p <- min(which(rspF >= maxF))
    # convert all times to ms from tenth of a ms
    #time from 10% to max force
    t10p <- min(which(rspF > 0.1 * maxF))
    #time from 50% to max force
    t50p <-
      t100p + min(which(rspF[t100p:length(rspF)] <= 0.5 * maxF))
    #maximum rate of change (from earlier F/s calcs)
    maxF1d <- max(rspF1d)
    #minimum rate of change
    minF1d <- min(rspF1d)
    #not sure what these two are
    t1dmax <- min(which(rspF1d >= maxF1d)) * 50
    t1dmin <- min(which(rspF1d <= minF1d)) * 50
    sumLine = c(
      sSpecies,
      snake,
      muscle,
      rater,
      MusMassg,
      i,
      meanF,
      maxF,
      t100p / 10,
      t10p / 10,
      t50p / 10,
      (t50p - t10p) / 10,
      maxF1d,
      minF1d,
      t1dmax / 10,
      t1dmin / 10,
      (t1dmin - t1dmax) / 10,
      sMAMU
    )
    #save these summary stats
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
}

#now do some formatting and add a little bit of information to the file
# we just wrote

df = read.csv("OutFiles/C4P/test/Couchii_C4P_Force.csv")
#df<-t(df[order(df$Pulse, df$Snake, df$Muscle),])

#remove outliers ----
library(ggplot2)
library(tidyverse)

# #get dataframe in format to actually plot (stats format =/= plot format)
# c4p_long <- df %>% 
#   pivot_longer(
#     cols = starts_with("X"),
#     names_to = "time",
#     values_to = "force"
#   ) %>% 
#   mutate(
#     time = as.numeric(gsub("X","", time))
#   )
# 
# #There's definitely a better way to do this, but I just manually split them so 
# # the colors were actually distinguishable and I can decide what the problem ones
# # are (i.e. which ones look like they tore and/or were normalized poorly)
# plot_IDs_1 <- ggplot(subset(c4p_long, Snake %in% c("CRF2630", "CRF2631", "CRF2633", 
#                                                    "CRF2669", "CRF2670", "CRF2671", 
#                                                    "CRF2672", "CRF2673", "CRF2674")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_1
# #2671 and 2669 don't seem to have really contracted
# 
# plot_IDs_2 <- ggplot(subset(c4p_long, Snake %in% c("CRF2676", "CRF2677", "CRF2678", 
#                                                    "CRF2679", "CRF2680", "CRF2681", 
#                                                    "CRF3051", "CRF3052", "CRF3055")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_2
# #all seem fine enough, some are slower to respond than others
# #actually 2677 should go based on how choppy it is; i think lots of these muscles were small and thus normalized fuzzily
# 
# plot_IDs_3 <- ggplot(subset(c4p_long, Snake %in% c("CRF3058", "CRF3059", "CRF3060", 
#                                                    "CRF3061", "CRF3064", "CRF3065", 
#                                                    "CRF3066", "CRF3069", "CRF3070")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_3
# #3066 way too big
# 
# plot_IDs_4 <- ggplot(subset(c4p_long, Snake %in% c("CRF3074", "CRF3211", "EJE164",  
#                                                    "EJE186", "CRF2675", "CRF3056", 
#                                                    "CRF3072")),
#                      aes(x = time, y = force, color = Snake)) +
#   geom_line() 
# plot_IDs_4
# #all pretty much fine
# 
# #filter out identified outliers and save file
# df <- df %>% 
#   filter(!Snake %in% c("CRF3066", "CRF2677", "CRF2671", "CRF2669"))
# 
# #write.csv(df,"OutFiles/C4P/test/Couchii_C4P_Force_Cleaned.csv")