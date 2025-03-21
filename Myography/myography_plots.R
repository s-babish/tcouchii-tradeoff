setwd("C:/Users/sdbab/OneDrive - University of Nevada, Reno/UNR/tcouchii/Myography")

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

#Plotting tetanus results -----
#First plot will plot all results separately ----
tet_force <- read.csv("Tcouchii_side_hustle/Couchii_tetanus/output/p3Tetanus-force.csv",
                      header = F)
#this file has the force output of the muscle, in N/g (grams of muscle), for 
# 1000 points pre-stimulus>1 to 30000 points after stimulus
# 10000 obs/sec, so 0.1 s before stimulus to 3 s after
# (also this file is huge so it takes quite a while)

#Normalized force = (measuredF - baseF)*9.80665/(Muscle Mass in grams)

colnames(tet_force) <- c("Species","Snake","Muscle","MAMU",paste0("t",1:30001))

tet_force_long <- tet_force %>% 
  pivot_longer(
    cols = starts_with("t"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("t","", time))
  )

#filter out outliers/invalid data based on waveform ----
#some muscles likely tore themselves or had other issues based on waveform, 
# and there's few enough runs I can filter those out manually instead of writing
# some sort of algorithm to recognize bad waveforms

#need to plot only 5 snakes at a time for colors to be distinguishable
#(i promise i'm aware of how incredibly cursed this method is, please forgive me)
plot_IDs_1 <- ggplot(subset(tet_force_long, 
                            Snake %in% c("CRF2630","CRF2631","CRF2633", "CRF2669","CRF2670")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_1
#lots of bad ones, 2669 and 2631 both look ripped and 2670 is a bit small

plot_IDs_2 <- ggplot(subset(tet_force_long, 
                            Snake %in% c("CRF2671", "CRF2672", "CRF2673", "CRF2674", "CRF2675")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_2
#all look fine, though the normalization is horrid
#2675 is the tall outlier, contrary to bobby's suggestion it's actually on the smaller side (11 mg)

plot_IDs_3 <- ggplot(subset(tet_force_long, 
                            Snake %in% c("CRF2676", "CRF2677", "CRF2678", "CRF2679", "CRF2680")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_3
#CRF2680 ripped

plot_IDs_4 <- ggplot(subset(tet_force_long, 
                            Snake %in% c("CRF2681", "CRF3051", "CRF3052", "CRF3055", "CRF3056")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_4
#all fine (except maybe 2681)

plot_IDs_5 <- ggplot(subset(tet_force_long, 
                            Snake %in% c("CRF3065", "CRF3066", "CRF3070", "CRF3074")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_5
#3074 totally wrong, 3065 and 3066 also look bad

plot_IDs_6 <- ggplot(subset(tet_force_long, 
                            Snake %in% c("CRF3058", "CRF3059", "CRF3060", "CRF3061", "CRF3064")),
                     aes(x = time, y = force, color = Snake)) +
  geom_line() 
plot_IDs_6
#3060 has weird decay
tet_force_long_sub <- tet_force_long %>% 
  filter(!Snake %in% c("CRF3060","CRF3074","CRF3065","CRF3066","CRF2680","CRF2669","CRF2631","CRF2670")) 
#%>% filter(Snake != "CRF2675")

#save this file
#write.csv(tet_force_long_sub, "data_processed/tetanus_no_outliers.csv")

#plot with MAMU color-coding -----
plot_MAMUs <- ggplot(tet_force_long_sub, aes(x = time, y = force, group = Snake, color = MAMU)) +
  geom_line() + 
  scale_color_viridis(option = "viridis")
plot_MAMUs

#Second plot will average results together ----
plot <- ggplot(tet_force_long, aes(x = time, y = force)) +
  geom_smooth()
plot

#alternative that fits an actual equation:
plot2 <- ggplot(tet_force_long,aes(x=time,y=force, group = Snake, color = MAMU))+
  geom_smooth(method="nls", 
              formula=y~1+Vmax*(1-exp(-x/tau)), # this is an nls argument
              method.args = list(start=c(tau=0.2,Vmax=2)), # this too
              se=FALSE) +
  scale_color_viridis(option = "viridis")
plot2

#Alternative alternative that fits bobby's equation:
#none of this works right now but i am not in the mood to figure it out
# plot3 <- ggplot(tet_force_long,aes(x=time,y=force, group = Snake, color = MAMU))+
#   geom_smooth(method="nls",
#               formula=y~A2 + (A1 * A2) / (1 + exp((log(x) - log(x0)) / dx)), # this is an nls argument
#               method.args = list(start=c(A1 = 1, A2 = 1, x0 = 400, dx = 1)), # this too
#               se=FALSE) +
#   scale_color_viridis(option = "viridis")
# plot3
# 
# drc_formula <- force~A2 + (A1 * A2) / (1 + exp((log(time) - log(x0)) / dx))

# tetanus_models <- tet_force_long %>%
#   group_by(Snake) %>%
#   do(fit = nlsLM(MMformula,., start = list(A2 = )))
# MMmodels %>% tidy(fit)


#**Variant including error bars ----

#Plotting tetanus derivative results ----
#basically the same process as the regular tetanus plots
tet_force1d <- read.csv("Tcouchii_side_hustle/Couchii_tetanus/output/p3Tetanus-force1d.csv",
header = F)

#  Low pass 200Hz filter -> pick every 30th entry & convert to per sec
#aka there are 600 F/s measurements

colnames(tet_force1d) <- c("Species","Snake","Muscle","MAMU",paste0("t",1:600))

tet_force1d_long <- tet_force1d %>% 
  filter(Snake != "CRF3074") %>% 
  pivot_longer(
    cols = starts_with("t"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("t","", time))
  )


plot <- ggplot(tet_force1d_long, aes(x = time, y = force, group = Snake)) +
  geom_line() + 
  scale_color_viridis(option = "viridis")
plot
#that one individual is so messed up - CRF3074; i'm excluding them
#this is another plot that will read way better once i can average it

