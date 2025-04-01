# C41P - Data analysis
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

setwd("C:/Users/sdbab/OneDrive - University of Nevada, Reno/UNR/tcouchii/Myography")
# Get Snake Info - Species, Genotype etc.
sinf = read.csv("Tcouchii_side_hustle/Scripts and datasheets/SnakeInfo-08.14.2019 (1).csv")

# Get Snake Muscle Masses
smm = read.csv("Tcouchii_side_hustle/Scripts and datasheets/SnakeSkeletalMuscleMasses-08.14.2019 (1).csv")
# reset all missing muscle mass values to -1.0smm
smm[is.na(smm)] <- 0.999

dname = "Tcouchii_side_hustle/Couchii_C4P"
# dname = "C:/Bobby/Data-CSV/p1c4pTest"
setwd(dname)

x = strsplit(dname, "/")

# Write / overwrite Column Headers for output summary file
# had to manually make this file
ofsum = paste("./output/", x[[1]][lengths(x)],  ".csv", sep = "")
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
# Output Force File
ofF <- "./output/c4p1-force.csv"
ofFhdr <- c(
  "Species",
  "MAMU",
  "MusMassg",
  "Pulse",
  "Snake",
  "Muscle",
  "Rater",
  (1:1500) / 10000
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
# Output Force first derivative File
#also had to manually make this one
ofF1d <- "./output/c4p1-force1d.csv"
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

# get all data files
files = list.files(path = ".", pattern = "csv")
q <- strsplit(files, "-")
# snakes <- unlist(lapply(q,'[[',1))
# muscles <- unlist(lapply(q,'[[',2))
# osum <- data.frame(files,snakes,muscles)

fC4P <- paste("./output/", "P1-C4PFiles",  ".csv", sep = "")
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

#a lot of the stuff in here references rows that don't exist in my snake info file which is at least part of the reason i'm having so many issues
for (file in files) {
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
  ofdata <- paste("./output/", snake, "-", muscle, ".csv", sep = "")
  raw = read.csv(file)
  # Select all rows with stimulus > 1
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
    # muscWeight = smm[as.factor(snake),muscle] / 1000
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
    maxF <- max(rspF)
    t100p <- min(which(rspF >= maxF))
    # convert all times to ms from tenth of a ms
    t10p <- min(which(rspF > 0.1 * maxF))
    t50p <-
      t100p + min(which(rspF[t100p:length(rspF)] <= 0.5 * maxF))
    maxF1d <- max(rspF1d)
    minF1d <- min(rspF1d)
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

df = read.csv("./output/c4p1-force.csv")
# soff <- off[order(off$Pulse,off$Snake.Muscle),]
# toff <- t(soff)
# write.csv(toff, "./output/c4p1-force-rpt.csv")
write.csv(t(df[order(df$Pulse, df$Snake, df$Muscle),]),
          "./output/c4p1-force-rpt.csv")
# Create reports by species and pulse
df <- df[order(df$Pulse, df$Species, df$Snake, df$Muscle), ]
for (p in unique(df$Pulse)) {
  for (s in levels(df$Species)) {
    write.csv(df[which(df$Pulse == p &
                          df$Species == s),], paste0("./output/p1C4P-Force-P",
                                                     p, "-", s, ".csv"))
  }
}

df = read.csv("./output/c4p1-force1d.csv")
write.csv(t(df[order(df$Pulse, df$Snake, df$Muscle),]),
          "./output/c4p1-force1d-rpt.csv")
df <- df[order(df$Pulse, df$Species, df$Snake, df$Muscle), ]
for (p in unique(df$Pulse)) {
  for (s in levels(df$Species[which(df$Pulse == p)])) {
    write.csv(df[which(df$Pulse == p &
                          df$Species == s),], paste0("./output/p1C4P-Force1d-P",
                                                     p, "-", s, ".csv"))
  }
}

# plot(rspF)
# plot(rspF1d,type="l")
# scatter.smooth(x=1:29, y=rspF1d)
