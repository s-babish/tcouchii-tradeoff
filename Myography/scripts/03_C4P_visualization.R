#libraries ----
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggforce)
library(ggpubr)
theme_set(theme_classic())

#load data ----
c4p_stats <- read.csv("OutFiles/C4P/Couchii_C4P_Metrics.csv")
c4p_raw <- read.csv("OutFiles/C4P/Couchii_C4P_Force_Cleaned.csv")

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
  ) %>% 
  drop_na()

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
  geom_smooth() +
  labs(x = "Time (s)", y = "Normalized Force")
plot_all

#Now recreate the trace plots for the first derivatives (Force/s) ----
#I suspect these plots will look better once I do my outlier removal
#Also I could calculate the derivative more frequently back in script 1
c4p_1d <- read.csv("OutFiles/C4P/Couchii_C4P_Force_1d_Cleaned.csv")

c4p_1d_long <- c4p_1d %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("X","", time))
  ) %>% 
  drop_na()

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

#read all in and code by genotype ----
#genotypes <- c("WT_elegans","WT_sirtalis","WT_hammondii","WT_atratus","LVNV","EPN","P","T")
genotypes <- c("LVNV","T","WT_sirtalis","WT_hammondii","EPN","WT_atratus")

#make storage matrix
# all_long <- matrix(ncol = 9)
# colnames(all_long) <- c("Species","Snake","Muscle","MAMU","MussMassg","Ind",
#                         "time","force","Genotype")
geno_means <- matrix(ncol = 4)
colnames(geno_means) <- c("avg","se","time","genotype")

for (type in genotypes) {
  foldername = type
  c4p_force <- read.csv(paste("OutFiles/C4P/",foldername,"/",foldername,
                              "_C4P_Force.csv",sep=""))
  
  colnames(c4p_force) <- c("Species","MAMU","MussMassg","Pulse","Snake","Muscle","Rater",
                           paste0("t",1:1500))
  c4p_force <- c4p_force %>% 
    mutate(
      Ind = paste(Snake,"_",Muscle,sep=""),
      Genotype = foldername
    ) %>% 
    filter(
      MussMassg != 0.000999,
      Snake != "CRF3110",
      Pulse == 1
      ) #i realize this goes nowhere for now, may use it later
  
  
  
  subset <- c4p_force[,8:1507] #only apply over the observations
  
  stats <- as.data.frame(t(sapply(subset, function(x) 
    c(avg = mean(x), se = sd(x)/sqrt(length(x)))))) %>% 
    mutate(
      time = as.numeric(gsub("t","",colnames(c4p_force)[8:1507]))/1000,
      genotype = foldername
    )
  
  geno_means <- rbind(geno_means,stats)
  
  
  #want to see what the individuals within the genotypes are doing
  c4p_long <- c4p_force %>% 
    pivot_longer(
      cols = starts_with("t"),
      names_to = "time",
      values_to = "force"
    ) %>% 
    mutate(
      time = as.numeric(gsub("t","", time)),
      time = time/1000
    )
  
  plot_in_geno <- ggplot(c4p_long, aes(x = time, y = force, group = Ind, color = Ind)) +
    geom_line() +
    ggtitle(foldername)
  print(plot_in_geno)
}

#write.csv(geno_means,"OutFiles/C4P/All_geno_avgs.csv",row.names=F)
plot_genotypes <- ggplot(geno_means[-1,], aes(x = time, y = avg, group = genotype, 
                                              fill = genotype, color = genotype)) +
  geom_line(lwd=1.5) #+ geom_ribbon(aes(ymin = avg - se, ymax = avg + se),alpha=0.25)
plot_genotypes
ggsave("Plots/transient_plots/geno_transient_comparisons.png",dpi=300)

#looking at metrics in box plot form ----
genotypes <- c("EPN","LVNV","T","WT_sirtalis","WT_atratus","WT_hammondii")
genotypes <- c("LVNV","T","WT_sirtalis","WT_hammondii")

c4p_all_long <- matrix(ncol = 4)
colnames(c4p_all_long) <- c("Snake","Genotype","variable","value")

for (type in genotypes) {
  foldername = type
  c4p_stats <- read.csv(paste("OutFiles/C4P/",foldername,"/",foldername,
                              "_C4P_Metrics.csv",sep=""))
c4p_long <- c4p_stats %>% 
  filter(
    MusMassg != 0.000999,
    Pulse == 1
  ) %>% 
  mutate(
    Genotype = foldername
  ) %>% 
  dplyr::select(Snake, BaseF.N.g.:DiffFChgMaxToMin.ms., Genotype) %>% 
  pivot_longer(cols = BaseF.N.g.:DiffFChgMaxToMin.ms., 
               names_to = "variable", values_to = "value")

  c4p_all_long <- rbind(c4p_all_long,c4p_long)
}
c4p_all_long <- c4p_all_long[-1,]
c4p_all_long$Genotype <- as.factor(c4p_all_long$Genotype)

for (page in 1:3) {
c4p_boxes <- ggplot(c4p_all_long, aes(Genotype, value, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, na.rm = T) +
  facet_wrap_paginate (. ~ variable, scales = 'free', shrink = T, nrow=2, ncol = 2, 
                       page = page) +
  xlab('') +
  ylab('') + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  geom_pwc(method = "dunn_test", p.adjust.method = "bonferroni", 
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("****", "***", "**", "*", "ns")),
           hide.ns = T, label = "p.adj.signif")
print(c4p_boxes)
ggsave(paste("Plots/transient_plots/Subset2_geno_metrics_",page,".png"),dpi=300)
}


#pairwise plots ----
#116BA8 for EPN color
#black for WT
#CC1E1D for LVNV
#1C8642 for T
#thin for line type to work
geno_means <- read.csv("OutFiles/C4P/All_geno_avgs.csv")
geno_means$genotype <- factor(geno_means$genotype, levels = c("WT_sirtalis","WT_elegans","T","LVNV","P","EPN"))

#geno_means <- geno_means[-1,]
png(filename="OutFiles/C4P/plots/WT_transient.png",width=3600,height=2700,res=500)
WT_trans <- ggplot(subset(geno_means, genotype %in% c("WT_sirtalis")), 
                        aes(x = time, y = avg, group = genotype,fill = genotype, 
                            color = genotype, linetype = genotype)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg - se, ymax = avg + se),alpha=0.25,color=NA) + 
  scale_linetype_manual(values = c("solid"),name = "Genotype", 
                        labels = c(expression(paste("WT ",italic("Th. sirtalis"))))) +
  scale_fill_manual(values = c("#000000"),
                    name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))))) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis")))),
                     values = c("#000000")) +
  labs(x = "Time (s)", y = "Scaled Force (N/g)",
       title = "Force of Transient Muscle Contractions") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16))

WT_trans
dev.off()

png(filename="OutFiles/C4P/plots/WT_T_transient.png",width=3600,height=2700,res=500)
WT_T_trans <- ggplot(subset(geno_means, genotype %in% c("WT_sirtalis","T")), 
                   aes(x = time, y = avg, group = genotype,fill = genotype, 
                       color = genotype, linetype = genotype)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg - se, ymax = avg + se),alpha=0.25,color=NA) + 
  scale_linetype_manual(values = c("solid","35"),name = "Genotype", 
                        labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                   expression(italic("Th. couchii")))) +
  scale_fill_manual(values = c("#000000","#1C8642"),
                    name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                                  expression(italic("Th. couchii")))) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                                   expression(italic("Th. couchii"))),
                     values = c("#000000","#1C8642")) +
  labs(x = "Time (s)", y = "Scaled Force (N/g)",
       title = "Force of Transient Muscle Contractions") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16))

WT_T_trans
dev.off()

png(filename="OutFiles/C4P/plots/WT_T_LVNV_transient.png",width=3600,height=2700,res=500)
WT_T_LVNV_trans <- ggplot(subset(geno_means, genotype %in% c("WT_sirtalis","T","LVNV")), 
                     aes(x = time, y = avg, group = genotype,fill = genotype, 
                         color = genotype, linetype = genotype)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg - se, ymax = avg + se),alpha=0.25, color = NA) + 
  scale_linetype_manual(values = c("solid","35","14"),name = "Genotype", 
                        labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                   expression(italic("Th. couchii")),
                                   expression(paste("LVNV ",italic("Th. sirtalis"))))) +
  scale_fill_manual(values = c("#000000","#1C8642","#CC1E1D"),
                    name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                                  expression(italic("Th. couchii")),
                                                  expression(paste("LVNV ",italic("Th. sirtalis"))))) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                                   expression(italic("Th. couchii")),
                                                   expression(paste("LVNV ",italic("Th. sirtalis")))),
                     values = c("#000000","#1C8642","#CC1E1D")) +
  labs(x = "Time (s)", y = "Scaled Force (N/g)",
       title = "Force of Transient Muscle Contractions") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16))

WT_T_LVNV_trans
dev.off()
