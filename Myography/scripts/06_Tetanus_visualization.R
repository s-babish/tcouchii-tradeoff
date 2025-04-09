#libraries ----
library(ggplot2)
library(tidyverse)
library(viridis)

#Plotting tetanus results -----
#First plot will plot all results separately ----
tet_force <- read.csv("OutFiles/Tetanus/test/Couchii_Tetanus_Force.csv")

#this file has the force output of the muscle, in N/g (grams of muscle), for 
# 1000 points pre-stimulus>1 to 30000 points after stimulus
# 10000 obs/sec, so 0.1 s before stimulus to 3 s after
# (also this file is huge so it takes quite a while)

#Normalized force = (measuredF - baseF)*9.80665/(Muscle Mass in grams)

colnames(tet_force) <- c("Species","Snake","Muscle","MAMU","MussMassg",
                         paste0("t",1:30000))

tet_force_long <- tet_force %>% 
  pivot_longer(
    cols = starts_with("t"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("t","", time))
  )



#plot with MAMU color-coding -----
plot_MAMUs <- ggplot(tet_force_long, aes(x = time, y = force, group = Snake, color = MAMU)) +
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
tet_force1d <- read.csv("OutFiles/Tetanus/test/Couchii_Tetanus_Force_1d.csv")

#It went through low pass 200Hz filter -> pick every 30th entry & convert to /s
#aka there are 600 F/s measurements

tet_force1d_long <- tet_force1d %>% 
  filter(Snake != "CRF3074") %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("X","", time))
  )


plot <- ggplot(tet_force1d_long, aes(x = time, y = force, group = Snake)) +
  geom_line() + 
  scale_color_viridis(option = "viridis")
plot
#that one individual is so messed up - CRF3074; i'm excluding them
#this is another plot that will read way better once i can average it

