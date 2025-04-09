#libraries ----
library(ggplot2)
library(tidyverse)
library(viridis)

#load data ----
c4p_stats <- read.csv("OutFiles/C4P/test/Couchii_C4P_Metrics.csv")
c4p_raw <- read.csv("OutFiles/C4P/test/Couchii_C4P_Force.csv")

#Comparing C4P stats between pulses (box plots) ----
#make long version of data:
c4p_long <- c4p_stats %>% 
  select(Snake, Pulse:DiffFChgMaxToMin.ms.) %>% 
  pivot_longer(cols = BaseF.N.g.:DiffFChgMaxToMin.ms., 
               names_to = "variable", values_to = "value")
c4p_long$Pulse <- as.factor(c4p_long$Pulse)

pulse_boxes <- ggplot(c4p_long, aes(Pulse, value, fill = Pulse)) +
  geom_boxplot(outlier.shape = NA, na.rm = T) +
  facet_wrap (. ~ variable, scales = 'free', shrink = T) +
  xlab('') +
  ylab('')
pulse_boxes
#all this tells me is that at least visually all the pulses are the same, so I won't waste time comparing 
# differences between pulses between MAMUs or IC50s

#Check for individual differences (mostly out of curiosity) ----

for (col in 7:17) {
  plot <- ggplot(c4p_stats, aes(x=Pulse, y=c4p_stats[,col], group=Snake, color=MAMU)) +
    geom_line(linewidth = 1) + 
    #geom_point() +
    ggtitle(paste0("Comparing inter-individual differences in ", colnames(c4p_stats)[col])) +
    ylab(paste0(colnames(c4p_stats)[col])) +
    theme_classic() +
    theme(axis.text = element_text(color="black"),panel.border = element_rect(fill = NA)) +
    scale_color_viridis(option = "viridis")
  print(plot)
}
#lmao i forgot how ugly these are (and that i probably won't ever use them)

#Plotting the force traces (pulses still separate) -----
#same normalized force as the tetanus results

#make dataframe long, change time value to make x axis actually work
c4p_long <- c4p_raw %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("X","", time))
  )

#Plot each pulse one at a time using subset()
pulse1 <- ggplot(subset(c4p_long, Pulse == 1), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Normalized Force")
pulse1

pulse2 <- ggplot(subset(c4p_long, Pulse == 2), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Normalized Force")
pulse2

pulse3 <- ggplot(subset(c4p_long, Pulse == 3), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Normalized Force")
pulse3

pulse4 <- ggplot(subset(c4p_long, Pulse == 4), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Normalized Force")
pulse4

#Combine traces and plot all 4 pulses together ----
#figure out what amount to add on to string all pulses together
end <- max(c4p_long$time)

#make all the time observation sequential so that they can easily be plotted together
#I tried to make this work with case_when and couldn't, so if anyone can that 
# would be super helpful (mostly for legibility but also bc I'm curious)
c4p_long <- c4p_long %>% 
  mutate(
    time_seq = ifelse(Pulse == 2, time + end, 
                      ifelse(Pulse == 3, time + 2*end,
                             ifelse(Pulse == 4, time + 3*end, time)))
  )

plot_all <- ggplot(c4p_long, aes(x = time_seq, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Normalized Force")
plot_all

#Now recreate the trace plots for the first derivatives (Force/s) ----
#I suspect these plots will look better once I do my outlier removal
#Also I could calculate the derivative more frequently back in script 1
c4p_1d <- read.csv("OutFiles/C4P/test/Couchii_C4P_Force_1d.csv")

c4p_1d_long <- c4p_1d %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("X","", time))
  )

#Plot each pulse one at a time using subset()
pulse1 <- ggplot(subset(c4p_1d_long, Pulse == 1), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Change in Force (delta F/s)")
pulse1

pulse2 <- ggplot(subset(c4p_1d_long, Pulse == 2), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Change in Force (delta F/s)")
pulse2

pulse3 <- ggplot(subset(c4p_1d_long, Pulse == 3), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Change in Force (delta F/s)")
pulse3

pulse4 <- ggplot(subset(c4p_1d_long, Pulse == 4), aes(x = time, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Change in Force (delta F/s)")
pulse4

#Combine traces and plot all 4 pulses together ----
#figure out what amount to add on to string all pulses together
end <- max(c4p_1d_long$time)

#make all the time observation sequential so that they can easily be plotted together
#I tried to make this work with case_when and couldn't, so if anyone can that 
# would be super helpful (mostly for legibility but also bc I'm curious)
c4p_1d_long <- c4p_1d_long %>% 
  mutate(
    time_seq = ifelse(Pulse == 2, time + end, 
                      ifelse(Pulse == 3, time + 2*end,
                             ifelse(Pulse == 4, time + 3*end, time)))
  )

plot_all <- ggplot(c4p_1d_long, aes(x = time_seq, y = force, group = Snake)) +
  geom_line() +
  labs(x = "Time (s)", y = "Change in Force (delta F/s)")
plot_all
