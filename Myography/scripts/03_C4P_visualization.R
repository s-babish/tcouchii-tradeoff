#libraries ----
library(ggplot2)
library(tidyverse)
library(viridis)

#Comparing C4P stats across pulses ----
#giving this a different name than it has in the other script to differentiate it
# from the other c4p file i'm gonna use
c4pstats <- read.csv("Tcouchii_side_hustle/Couchii_C4P/output/Couchii_C4P.csv")

for (col in 7:17) {
    plot <- ggplot(c4pstats, aes(x=Pulse, y=c4pstats[,col], group=Snake, color=MAMU)) +
    geom_line(linewidth = 1) + 
    #geom_point() +
    ggtitle(paste0("Comparing inter-individual differences in ", colnames(c4pstats)[col])) +
    ylab(paste0(colnames(c4pstats)[col])) +
      theme_classic() +
      theme(axis.text = element_text(color="black"),panel.border = element_rect(fill = NA)) +
      scale_color_viridis(option = "viridis")
    print(plot)
}

#Plotting C4P (transient force) results, pulses plotted individually -----
c4p_raw <- read.csv("Tcouchii_side_hustle/Couchii_C4P/output/c4p1-force.csv")

#same normalized force as the tetanus results

c4p_long <- c4p_raw %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("X","", time))
  )
#remove weird ones (based on pulse 1 waveform plots)
plot_IDs_1 <- ggplot(subset(c4p_long, Snake %in% c("CRF2630", "CRF2631", "CRF2633", "CRF2669", "CRF2670", "CRF2671", "CRF2672", "CRF2673", "CRF2674")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_1
#2671 and 2669 don't seem to have really contracted

plot_IDs_2 <- ggplot(subset(c4p_long, Snake %in% c("CRF2676", "CRF2677", "CRF2678", "CRF2679", "CRF2680", "CRF2681", "CRF3051", "CRF3052", "CRF3055")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_2
#all seem fine enough, some are slower to respond than others
#actually 2677 should go based on how choppy it is; i think lots of these muscles were small and thus normalized fuzzily

plot_IDs_3 <- ggplot(subset(c4p_long, Snake %in% c("CRF3058", "CRF3059", "CRF3060", "CRF3061", "CRF3064", "CRF3065", "CRF3066", "CRF3069", "CRF3070")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_3
#3066 way too big

plot_IDs_4 <- ggplot(subset(c4p_long, Snake %in% c("CRF3074", "CRF3211", "EJE164",  "EJE186", "CRF2675", "CRF3056", "CRF3072")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_4
#all pretty much fine

#remove outliers from visual examination
c4p_long_sub <- c4p_long %>% 
  filter(!Snake %in% c("CRF3066", "CRF2677", "CRF2671", "CRF2669"))

#save file
#write.csv(c4p_long_sub, "data_processed/c4p_no_outliers.csv")

c4p1 <- c4p_long_sub %>% 
  filter(Pulse == 1) 

plot1 <- ggplot(c4p1, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot1

c4p2 <- c4p_long_sub %>% 
  filter(Pulse == 2)

plot2 <- ggplot(c4p2, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot2

c4p3 <- c4p_long_sub %>% 
  filter(Pulse == 3)

plot3 <- ggplot(c4p3, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot3

c4p4 <- c4p_long_sub %>% 
  filter(Pulse == 4)

plot4 <- ggplot(c4p4, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot4