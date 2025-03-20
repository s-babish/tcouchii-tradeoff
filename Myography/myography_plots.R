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

c4p1 <- c4p_long %>% 
  filter(Pulse == 1) %>% 
  filter(Snake!= "CRF3600") %>% 
  filter(Snake != "CRF3066")
c4p2 <- c4p_long %>% 
  filter(Pulse == 2)
c4p3 <- c4p_long %>% 
  filter(Pulse == 3)
c4p4 <- c4p_long %>% 
  filter(Pulse == 4)

plot1 <- ggplot(c4p1, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot1
#another weird one again; also notable that the force is a normal amount in this
# file, i need to look into why
c4p1[c4p1$force > 8, ]
#CRF3600 and 3066, so of course multiple different ones this time

c4p2 <- c4p_long %>% 
  filter(Pulse == 2)%>% 
  filter(Snake!= "CRF3600") %>% 
  filter(Snake != "CRF3066")

plot2 <- ggplot(c4p2, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot2

c4p3 <- c4p_long %>% 
  filter(Pulse == 3)
  filter(Snake!= "CRF3600") %>% 
  filter(Snake != "CRF3066")

plot3 <- ggplot(c4p3, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot3

pc4p4 <- c4p_long %>% 
  filter(Pulse == 4)%>% 
  filter(Snake!= "CRF3600") %>% 
  filter(Snake != "CRF3066")

plot4 <- ggplot(c4p4, aes(x = time, y = force, group = Snake)) +
  geom_line() 
plot4

#I think bobby might have averaged this stuff together? unsure
#also these plots won't look good until i figure out the averaging thing

#Plotting C4P (transient force) results, pulses plotted sequentially ----

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


plot <- ggplot(tet_force_long, aes(x = time, y = force, group = Snake, color = MAMU)) +
  geom_line() + 
  scale_color_viridis(option = "viridis")
plot

#figure out who the really strong individual was
tet_force_long[tet_force_long$force > 10000, ]
#CRF2675

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
plot3 <- ggplot(tet_force_long,aes(x=time,y=force, group = Snake, color = MAMU))+
  geom_smooth(method="nls", 
              formula=A2 + (A1 * A2) / (1 + exp((log(x) - log(x0)) / dx)), # this is an nls argument
              method.args = list(start=c(A1 = 1, A2 = 1, x0 = 1, dx = 1)), # this too
              se=FALSE) + 
  scale_color_viridis(option = "viridis")
plot3


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
#Rheobase plots (separated by pulse length) -----
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
