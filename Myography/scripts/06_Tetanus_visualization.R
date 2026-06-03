#libraries ----
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)
library(ggforce)
theme_set(theme_classic())

#set folder name
foldername <- "WT_atratus"

#Plotting tetanus results -----
#First plot will plot all results separately ----
tet_force <- read.csv(paste("OutFiles/Tetanus/",foldername,"/",foldername,"_Tetanus_Force.csv",sep=""))

#this file has the force output of the muscle, in N/g (grams of muscle), for 
# 1000 points pre-stimulus>1 to 30000 points after stimulus
# 10000 obs/sec, so 0.1 s before stimulus to 3 s after
# (also this file is huge so it takes quite a while)
#Normalized force = (measuredF - baseF)*0.00980665/(Muscle Mass in grams)

colnames(tet_force) <- c("Species","Snake","Muscle","MAMU","MussMassg",
                         paste0("t",1:30000))

tet_force_long <- tet_force %>% 
  mutate(
    Ind = paste(Snake,"_",Muscle,sep="")
  ) %>% 
  pivot_longer(
    cols = starts_with("t"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("t","", time))
  )



#plot with MAMU color-coding -----
plot_MAMUs <- ggplot(tet_force_long, aes(x = time, y = force, group = Ind, color = MAMU)) +
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
#this might work better if i scale the tetanus force outputs
plot3 <- ggplot(tet_force_long,aes(x=time,y=force, group = Snake, color = MAMU))+
  geom_smooth(method="nls",
              formula=y~A2 + (A1 - A2) / (1 + exp((log(x) - log(exp(p0))) / dx)), # this is an nls argument
              method.args = list(start=c(A1 = 0, A2 = 3000, p0 = 6, dx = 1)), # this too
              data = tet_force_long,
              se=FALSE) +
  scale_color_viridis(option = "viridis") +
  scale_x_continuous(limits = c(0,10000))
plot3

#still isn't fitting the ones further off from the given start values but it's
# an improvement for sure (and limiting x-vals helped a lot)

#**Variant including error bars ----

#Plotting tetanus derivative results ----
#basically the same process as the regular tetanus plots
tet_force1d <- read.csv(paste("OutFiles/Tetanus/",foldername,"/",foldername,"_Tetanus_Force_1d.csv",sep=""))

#It went through low pass 200Hz filter -> pick every 30th entry & convert to /s
#aka there are 600 F/s measurements

tet_force1d_long <- tet_force1d %>% 
  mutate(
    Ind = paste(Snake,"_",Muscle,sep="")
  ) %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "time",
    values_to = "force"
  ) %>% 
  mutate(
    time = as.numeric(gsub("X","", time))
  )


plot <- ggplot(tet_force1d_long, aes(x = time, y = force, group = Ind)) +
  geom_line() + 
  scale_color_viridis(option = "viridis")
plot
#that one individual is so messed up - CRF3074; i'm excluding them
#this is another plot that will read way better once i can average it

#testing metrics ----
tet_metrics <- read.csv(paste("OutFiles/Tetanus/",foldername,"/",foldername,"_Tetanus_Metrics.csv",sep=""))
colnames(tet_metrics)
tet_metrics$Snake
tet_force$Snake

for (row in 1:nrow(tet_force)) {
  plot(1:30000,tet_force[row,c(6:ncol(tet_force))])
  abline(h=tet_metrics[row,5])
  abline(h=tet_metrics[row,14])
  title(main = tet_force$Ind[row])
}

#read all in and code by genotype ----
genotypes <- c("WT_elegans","WT_sirtalis","WT_hammondii","WT_atratus","LVNV","EPN","P","T")
#genotypes <- c("WT_sirtalis","LVNV","T")

geno_means <- matrix(ncol = 4)
colnames(geno_means) <- c("avg","se","time","genotype")

for (type in genotypes) {
  foldername = type
  tet_force <- read.csv(paste("OutFiles/Tetanus/",foldername,"/",foldername,
                              "_Tetanus_Force.csv",sep=""))
  
  colnames(tet_force) <- c("Species","Snake","Muscle","MAMU","MussMassg",
                           paste0("t",1:30000))
  tet_force <- tet_force %>% 
    mutate(
    Ind = paste(Snake,"_",Muscle,sep=""),
    Genotype = foldername
  ) %>% 
    filter(MussMassg != 0.000999)
  
  subtet <- tet_force[,6:30005]
  
  stats <- (sapply(subtet, function(x) 
    c(avg = mean(x), se = sd(x)/sqrt(length(x)))))
    
    stats <- as.data.frame(t(stats)) %>% 
        mutate(
        time = seq(0,2.9999,by=0.0001),
        genotype = foldername
      )
  
  geno_means <- rbind(geno_means,stats)

}

#write.csv(geno_means,"OutFiles/Tetanus/sub_geno_avgs.csv",row.names=F) #just 3 from talk
write.csv(geno_means, "OutFiles/Tetanus/All_geno_avgs.csv",row.names=F)

plot_genotypes <- ggplot(geno_means[-1,], aes(x = time, y = avg, group = genotype,
                                              fill = genotype, color = genotype)) +
  geom_line() 
#+
 # geom_ribbon(aes(ymin = avg - se, ymax = avg + se),alpha=0.5)
plot_genotypes
ggsave("OutFiles/Tetanus/geno_tetanus_comparisons.png",dpi=300)

#looking at metrics in box plot form ----
tet_all_long <- matrix(ncol = 4)
colnames(tet_all_long) <- c("Snake","Genotype","variable","value")

for (type in genotypes) {
  foldername = type
  tet_stats <- read.csv(paste("OutFiles/Tetanus/",foldername,"/",foldername,
                              "_Tetanus_Metrics.csv",sep=""))
  tet_long <- tet_stats %>%
    filter(
      MusMassg != 0.000999
    ) %>%
    mutate(
      Genotype = foldername
    ) %>%
    dplyr::select(Snake, BaseF.N.g.:DiffFMaxToEnd, Genotype) %>%
    pivot_longer(cols = BaseF.N.g.:DiffFMaxToEnd,
                 names_to = "variable", values_to = "value")

  tet_all_long <- rbind(tet_all_long,tet_long)
}
tet_all_long <- tet_all_long[-1,]
tet_all_long$Genotype <- as.factor(tet_all_long$Genotype)

for (page in 1:6) {
  tet_boxes <- ggplot(tet_all_long, aes(Genotype, value, fill = Genotype)) +
    geom_boxplot(outlier.shape = NA, na.rm = T) +
    facet_wrap_paginate (. ~ variable, scales = 'free', shrink = T, nrow=1, ncol = 2,
                         page = page) +
    xlab('') +
    ylab('') +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    geom_pwc(method = "dunn_test", p.adjust.method = "bonferroni",
             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                symbols = c("****", "***", "**", "*", "ns")),
             hide.ns = T, label = "p.adj.signif")
  print(tet_boxes)
  ggsave(paste("OutFiles/Tetanus/plots/geno_metrics_",page,".png"),dpi=300)
}

#twitch-tet comps ----
twitch_tet_comparisons <- read.csv("OutFiles/Tetanus/twitch_tet_comps.csv")
tet_boxes <- ggplot(twitch_tet_comparisons, aes(Genotype, twitch_tet_ratio, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, na.rm = T) +
  xlab('') +
  ylab('') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  geom_pwc(method = "dunn_test", p.adjust.method = "bonferroni",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                              symbols = c("****", "***", "**", "*", "ns")),
           hide.ns = T, label = "p.adj.signif")
print(tet_boxes)
ggsave(paste("OutFiles/Tetanus/plots/geno_metrics_",page+1,".png"),dpi=300)

#pairwise plots ----
df <- read.csv("OutFiles/Tetanus/sub_geno_avgs.csv")
#thin data so line type is legible
df.new = df[seq(1, nrow(df), 5), ]
df.new$genotype <- factor(df.new$genotype, levels = c("WT_sirtalis","T","LVNV"))

png(filename="OutFiles/Tetanus/plots/WT_tetanic.png",width=4000,height=2700,res=500)
WT_tet <- ggplot(subset(df.new, genotype %in% c("WT_sirtalis")), 
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
       title = "Force of Tetanic Muscle Contractions") +
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

WT_tet
dev.off()

png(filename="OutFiles/Tetanus/plots/WT_T_tetanic.png",width=4000,height=2700,res=500)
WT_T_tet <- ggplot(subset(df.new, genotype %in% c("WT_sirtalis","T")), 
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
       title = "Force of Tetanic Muscle Contractions") +
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

WT_T_tet
dev.off()

png(filename="OutFiles/Tetanus/plots/WT_T_LVNV_tetanic.png",width=4000,height=2700,res=500)
WT_T_LVNV_tet <- ggplot(subset(df.new, genotype %in% c("WT_sirtalis","T","LVNV")), 
                   aes(x = time, y = avg, group = genotype,fill = genotype, 
                       color = genotype,linetype = genotype)) +
  geom_line() + 
  geom_ribbon(aes(ymin = avg - se, ymax = avg + se),alpha=0.25,color=NA) + 
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
       title = "Force of Tetanic Muscle Contractions") +
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

WT_T_LVNV_tet
dev.off()