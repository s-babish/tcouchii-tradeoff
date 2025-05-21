library(tidyverse)
library(ggplot2)
library(viridis)
library(ggforce)

#Scaled force vs time (separated by pulse length) -----
#vector with all pulse length options for iterating (excluding files with no or few snakes)
pulse_l <- c(50,100,200,500,1000,10000,50000)
i=1

keydf <- read.csv("OutFiles/Rheobase/Sigmoidal/test/Sigmoidal-rpt.csv")
keydf <- keydf[,c(7,11)]

for (i in 1:length(pulse_l)) {
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
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
      Snake = sub("_.*", "", Snake),
      MAMU = 0
    ) %>% 
    rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
  
  plot <- ggplot(df_long, aes(x = Dose, y = scaled_force, group = Snake, color = MAMU)) +
    geom_smooth(se = F) + 
    #geom_point() +
    #geom_line(aes(group = Snake)) +
    scale_color_viridis(option = "viridis") +
    scale_x_continuous(labels = dose_labels, breaks = 0:15) +
    ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
    xlab("Current (mA)") + 
    ylab("Proportion of maximal force")
  print(plot)
  
}

#actually fit the curve to the data ----
i=1
for (i in 1:length(pulse_l)) {
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile)
  df <- df[,c(-1)] #remove index line
  dose_labels = df$Dose
  df_long <- df %>%
    pivot_longer(
      cols = -Dose,
      names_to = "Snake",
      values_to = "scaled_force"
    ) %>% 
    mutate(
      Snake = sub("_.*", "", Snake),
      MAMU = 0
    ) %>% 
    rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
  
  plot <- ggplot(df_long, aes(x = Dose, y = scaled_force, group = Snake, color = MAMU)) + 
    geom_point() + 
    geom_smooth(method = "nls", 
                method.args = list(formula = y ~ (-1)/(1 + exp((log(x)-log(x0))/dx)) + 1,
                                   start = list(x0 = 60 , dx = 1)
                                    # ,
                                    # control = nls.lm.control(maxiter = 1024, maxfev = 1024),
                                    # lower = c(0,-Inf), upper = c(600,Inf), algorithm = "port"
                                   ), 
                data = df_long,
                se = FALSE) +
    scale_color_viridis(option = "viridis") +
    ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
    xlab("Current (mA)") + 
    ylab("Proportion of maximal force")
  print(plot)
  
}


#force/mA vs pulse width*mA ----
#note that I don't think this is really a proper plot bc these are the scaled
# force values, but I do not want to have to deal with modifying the python script 
# just for these plots at this point, so this is what we have for now

for (j in 1:length(pulse_l)) {
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[j],".csv", sep=""))
  df <- read.csv(infile)
  df <- df[,c(-1)] #remove index line
  #multiple force in each column by mA
  force_df <- mapply('/',df[3:ncol(df)],df[2])
  force_df[1,] <- 0
  force_df_long <- as.data.frame(force_df) %>% 
    mutate(
      mA_s = df$Dose * pulse_l[j]
    ) %>% 
    pivot_longer(
      cols = -mA_s,
      names_to = "Snake",
      values_to = "g_mA"
    ) %>% 
    mutate(
      Snake = sub("_.*", "", Snake),
      MAMU = 0
    ) %>% 
    rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
  
  plot <- ggplot(force_df_long, aes(x = mA_s, y = g_mA, group = Snake, color = MAMU)) +
    #geom_point() +
    geom_line(aes(group = Snake)) +
    scale_color_viridis(option = "viridis") +
    scale_x_continuous(labels = dose_labels, breaks = 0:15) +
    ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
    xlab("Total current (mA*s)") + 
    ylab("Proportion of maximal force per mA")
  print(plot)
  
}
#I suspect these plots will make more sense if/when I do it with the not-scaled
# force data, I think the sub-0 values are messing with it


#pulse width vs x0 ------
#renaming to make more sense for what we're doing here
rheobase_sigmoidal_plot <- read.csv("OutFiles/Rheobase/Sigmoidal/test/Sigmoidal-rpt.csv")
#remove rows we don't want
rheobase_sigmoidal_plot <- rheobase_sigmoidal_plot[-1,-1]

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

#ANOVA and Tukey's confirm the pulses have statistically equivalent x0 values

#group individuals
#cursed transformation to show all pulse lengths equidistant
rheobase_sigmoidal_plot$pulse_length <- as.numeric(rheobase_sigmoidal_plot$pulse_length)

rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=x0, color=X),
                                    show.legend=F) +
  geom_point(show.legend = F) +
  scale_x_continuous(labels = pulse_l, breaks = 1:7) +
  geom_line(show.legend=F)

# Integrate MAMU into plots ----
rheobase_sigmoidal_plot %>% 
  ggplot(aes(MAMU, x0, color = as.factor(pulse_length)))+
  geom_point()

rheobase_sigmoidal_plot %>% 
  ggplot(aes(MAMU, x0))+
  geom_point() + 
  facet_wrap(~as.factor(pulse_length))

#IC50 plots (data formatted in analyses script for now) -----
plot(x0 ~ IC50, data = IC50dat)

#plots for norm and bobby ----
#storage dataframe
all_rheo <- data.frame(
  Snake = character(),
  scaled_force = double(),
  pulse_length = integer(),
  MAMU = double()
)

for (i in 1:length(pulse_l)) {
  infile <- toString(paste("OutFiles/Rheobase/Force_scaled/TPA_",pulse_l[i],".csv", sep=""))
  df <- read.csv(infile)
  df <- df[,c(-1)] #remove index line
  dose_labels = df$Dose
  df_long <- df %>%
    pivot_longer(
      cols = -Dose,
      names_to = "Snake",
      values_to = "scaled_force"
    ) %>% 
    mutate(
      Snake = sub("_.*", "", Snake),
      pulse_length = pulse_l[i],
      MAMU = 0
    ) %>% 
    rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") %>% #pull in MAMU
    mutate(
      pulse_length = as.factor(pulse_length)
    )
  all_rheo <- rbind(all_rheo,df_long)
}

for (page in 1:6) {
plot <- ggplot(all_rheo, aes(x = Dose, y = scaled_force, group = pulse_length, color = pulse_length)) +
  geom_point() + 
  facet_wrap_paginate(~ Snake, ncol = 3, nrow = 2, page = page) +
  geom_smooth(method = "nls",
              method.args = list(formula = y ~ (-1)/(1 + exp((log(x)-log(x0))/dx)) + 1,
                                 start = list(x0 = 60 , dx = 1)
                                 # ,
                                 # control = nls.lm.control(maxiter = 1024, maxfev = 1024),
                                 # lower = c(0,-Inf), upper = c(600,Inf), algorithm = "port"
              ),
              data = all_rheo,
              se = FALSE) +
  scale_color_viridis(discrete = T) +
  xlab("Current (mA)") + 
  ylab("Proportion of maximal force")
print(plot)
}
