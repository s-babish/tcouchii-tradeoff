library(tidyverse)
library(ggplot2)
library(viridis)

#Scaled force vs time (separated by pulse length) -----
#vector with all pulse length options for iterating (excluding files with no or few snakes)
pulse_l <- c(50,100,200,500,1000,10000,50000)
i=1

keydf <- read.csv("TPA-sigmoidal-rpt.csv")
keydf <- keydf[,c(7,11)]

for (i in 1:length(pulse_l)) {
  infile <- toString(paste("TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile)
  df <- df[,c(-1)] #remove index line
  dose_labels = df$Dose
  df_long <- df %>%
    mutate (
      Dose = 0:15 #this is to make the plot look better bc ggplot is fighting me
    ) %>% 
    pivot_longer(
      cols = -Dose,
      names_to = "Snake",
      values_to = "scaled_force"
    ) %>% 
    mutate(
      Snake = sub("_[^_]+$", "", Snake),
      MAMU = 0
    ) %>% 
    rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
  
  plot <- ggplot(df_long, aes(x = Dose, y = scaled_force, group = Snake, color = MAMU)) +
    geom_smooth(se = F) + 
    scale_color_viridis(option = "viridis") +
    scale_x_continuous(labels = dose_labels, breaks = 0:15) +
    ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
    xlab("Current (mA)") + 
    ylab("Proportion of maximal force")
  print(plot)
  
}

#pulse width vs x0 ------
#(data formatted in rheobase_analyses.R)
rheobase_sigmoidal_plot <- rheobase_sigmoidal

#general view
rheobase_sigmoidal_plot$pulse_length <- as.factor(rheobase_sigmoidal_plot$pulse_length)

rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=x0, fill=pulse_length)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Rheobase x0 vs pulse width") +
  xlab("Pulse width (us)")

#group individuals
#cursed transformation to show all pulse lengths equidistant
rheobase_sigmoidal_plot$pulse_length <- as.numeric(rheobase_sigmoidal_plot$pulse_length)

rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=x0, color=X),
                                    show.legend=F) +
  geom_point(show.legend = F) +
  scale_x_continuous(labels = pulse_l, breaks = 1:7) +
  geom_line(show.legend=F)

