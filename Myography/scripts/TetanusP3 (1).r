# TetanusP3 - Data analysis
#   Protocol 3: One 2sec pulse
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
setwd("C:/Users/sdbab/OneDrive - University of Nevada, Reno/UNR/tcouchii/Myography")
# Get Snake Info - Species, Genotype etc.
sinf = read.csv("Tcouchii_side_hustle/Scripts and datasheets/SnakeInfo-08.14.2019 (1).csv")

# Get Snake Muscle Masses
smm = read.csv("Tcouchii_side_hustle/Scripts and datasheets/SnakeSkeletalMuscleMasses-08.14.2019 (1).csv")
# reset all missing muscle mass values to 0.999
smm[is.na(smm)] <- 0.999

dname = "Tcouchii_side_hustle/Couchii_tetanus"
# dname = "C:/Bobby/Data-CSV/p3Tetanus"
setwd(dname)

x = strsplit(dname, "/")

# Write / overwrite Column Headers for output summary file
ofsum = paste("./output/", x[[1]][lengths(x)],  ".csv", sep = "")
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
#
ofF <- "./output/p3Tetanus-force.csv"
ofFhdr <-
  c(
    "Species",
    "Snake",
    "Muscle",
    "MAMU",
    "MusMassg",
    (1:30000) / 10000
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
#
ofF1d <- "./output/p3Tetanus-force1d.csv"
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
#
# get all data files
files = list.files(path = ".", pattern = "csv")
q <- strsplit(files, "-")
# snakes <- unlist(lapply(q,'[[',1))
# muscles <- unlist(lapply(q,'[[',2))
# osum <- data.frame(files,snakes,muscles)

fTet <- paste("./output/", "P3-TetanusFiles",  ".csv", sep = "")
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

for (file in files) {
  y <- strsplit(file, "-")
  snake <- trimws(y[[1]][1])
  # get Snake Info
  sSpecies <- as.character(sinf$Species[which(sinf$Snake == snake)])
  sMAMU <- as.character(sinf$MAMU[which(sinf$Snake == snake)])
  # get Muscle Mass
  muscle <- trimws(y[[1]][2])
  sMusMassg <- smm[which(smm$SnakeID == snake), muscle] / 1000
  # write.table(
  #   t(c(file, snake, muscle, sMusMassg)),
  #   file = fTet,
  #   append = TRUE,
  #   quote = TRUE,
  #   sep = ",",
  #   eol = "\n",
  #   na = "NA",
  #   dec = ".",
  #   row.names = FALSE,
  #   col.names = FALSE
  # )
  
  ofdata <- paste("./output/", snake, "-", muscle, ".csv", sep = "")
  raw = read.csv(file)
  # Select all rows with stimulus > 1
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
  # write.table(
  #   t(oFline),
  #   file = ofF,
  #   append = TRUE,
  #   quote = TRUE,
  #   sep = ",",
  #   eol = "\n",
  #   na = "NA",
  #   dec = ".",
  #   row.names = FALSE,
  #   col.names = FALSE
  # )
  #  Low pass 200Hz filter -> pick every 30th entry & convert to per sec
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


#this splits by species which is unnecessary for me obviously

# df = read.csv("./output/p3tetanus-force.csv")
# df <- df[order(df$Species, df$Snake, df$Muscle),]
# for (s in unique(df$Species)) {
#   write.csv(df[which(df$Species == s), ], paste0("./output/p3Tetanus-Force-", s, ".csv"))
# }
# 
# df = read.csv("./output/p3tetanus-force1d.csv")
# df <- df[order(df$Species, df$Snake, df$Muscle),]
# for (s in unique(df$Species)) {
#   write.csv(df[which(df$Species == s), ], paste0("./output/p3Tetanus-Force1d-", s, ".csv"))
# }

#
# plot(rspF)
# plot(rspF1d,type="l")
# scatter.smooth(x=1:29, y=rspF1d)
