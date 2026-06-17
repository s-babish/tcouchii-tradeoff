library(tidyverse)
library(ggplot2)
library(viridis)
library(ggforce)
theme_set(theme_classic())

#set according to dataset 
foldername <- "WT_elegans"
#Scaled force vs time (separated by pulse length) -----
#vector with all pulse length options for iterating (excluding files with no or few snakes)
pulse_l <- c(50,100,200,500,1000,10000,50000)

keydf <- read.csv(toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv", sep="")))
keydf <- keydf[,c(10,14)]

# for (i in 1:length(pulse_l)) {
#   infile <- toString(paste("OutFiles/Rheobase/",foldername,"/Force_scaled/Rheo_",pulse_l[i],".csv", sep=""))
#   df <- read.csv(infile)
#   df <- df[-17,c(-1)] #remove index line
#   dose_labels = df$Dose
#   df_long <- df %>%
#     mutate (
#       Dose = 0:15 #this is to make the plot look better bc ggplot is fighting me
#     ) %>% 
#     pivot_longer(
#       cols = -Dose,
#       names_to = "Snake",
#       values_to = "scaled_force"
#     ) %>% 
#     mutate(
#       Snake = sub("_.*", "", Snake),
#       MAMU = 0
#     ) %>% 
#     rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
#   
#   plot <- ggplot(df_long, aes(x = Dose, y = scaled_force, group = Snake, color = MAMU)) +
#     #geom_smooth(se = F) + 
#     geom_point() +
#     geom_line(aes(group = Snake)) +
#     scale_color_viridis(option = "viridis") +
#     scale_x_continuous(labels = dose_labels, breaks = 0:15) +
#     ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
#     xlab("Current (mA)") + 
#     ylab("Proportion of maximal force")
#   print(plot)
#   
# }
# 
# #actually fit the curve to the data ----
# for (i in 1:length(pulse_l)) {
#   infile <- toString(paste("OutFiles/Rheobase/",foldername,"/Force_scaled/Rheo_",pulse_l[i],".csv", sep=""))
#   df <- read.csv(infile)
#   df <- df[,c(-1)] #remove index line
#   dose_labels = df$Dose
#   df_long <- df %>%
#     pivot_longer(
#       cols = -Dose,
#       names_to = "Snake",
#       values_to = "scaled_force"
#     ) %>% 
#     mutate(
#       Snake = sub("_.*", "", Snake),
#       MAMU = 0
#     ) %>% 
#     rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
#   
#   plot <- ggplot(df_long, aes(x = Dose, y = scaled_force, group = Snake, color = MAMU)) + 
#     geom_point() + 
#     geom_smooth(method = "nls", 
#                 method.args = list(formula = y ~ (-1)/(1 + exp((log(x)-log(x0))/dx)) + 1,
#                                    start = list(x0 = 60 , dx = 1)
#                                     # ,
#                                     # control = nls.lm.control(maxiter = 1024, maxfev = 1024),
#                                     # lower = c(0,-Inf), upper = c(600,Inf), algorithm = "port"
#                                    ), 
#                 data = df_long,
#                 se = FALSE) +
#     scale_color_viridis(option = "viridis") +
#     ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
#     xlab("Current (mA)") + 
#     ylab("Proportion of maximal force")
#   print(plot)
#   
# }
# 
# 
# #force/mA vs pulse width*mA ----
# #note that I don't think this is really a proper plot bc these are the scaled
# # force values, but I do not want to have to deal with modifying the python script 
# # just for these plots at this point, so this is what we have for now
# 
# for (j in 1:length(pulse_l)) {
#   infile <- toString(paste("OutFiles/Rheobase/",foldername,"/Force_scaled/Rheo_",pulse_l[j],".csv", sep=""))
#   df <- read.csv(infile)
#   df <- df[,c(-1)] #remove index line
#   #multiple force in each column by mA
#   force_df <- mapply('/',df[3:ncol(df)],df[2])
#   force_df[1,] <- 0
#   force_df_long <- as.data.frame(force_df) %>% 
#     mutate(
#       mA_s = df$Dose * pulse_l[j]
#     ) %>% 
#     pivot_longer(
#       cols = -mA_s,
#       names_to = "Snake",
#       values_to = "g_mA"
#     ) %>% 
#     mutate(
#       Snake = sub("_.*", "", Snake),
#       MAMU = 0
#     ) %>% 
#     rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") #pull in MAMU
#   
#   plot <- ggplot(force_df_long, aes(x = mA_s, y = g_mA, group = Snake, color = MAMU)) +
#     #geom_point() +
#     geom_line(aes(group = Snake)) +
#     scale_color_viridis(option = "viridis") +
#     scale_x_continuous(labels = dose_labels, breaks = 0:15) +
#     ggtitle(toString(paste("Rheobase, pulse length ", pulse_l[i], "us"))) +
#     xlab("Total current (mA*s)") + 
#     ylab("Proportion of maximal force per mA")
#   print(plot)
#   
# }
# #I suspect these plots will make more sense if/when I do it with the not-scaled
# # force data, I think the sub-0 values are messing with it


#pulse width vs x0/EC50 ------
#renaming to make more sense for what we're doing here
rheobase_sigmoidal_plot <- read.csv(toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv",sep="")))
#remove rows we don't want
rheobase_sigmoidal_plot <- rheobase_sigmoidal_plot[-1,-1]

#general view
rheobase_sigmoidal_plot$pulse_length <- as.factor(rheobase_sigmoidal_plot$pulse_length)

ec50 <- rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=x0, fill=pulse_length)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste(foldername,"Rheobase EC50 vs pulse width",sep=" ")) +
  xlab("Pulse width (us)")

ggsave(paste("OutFiles/Rheobase/",foldername,"/plots/",foldername,"_EC50.png",sep=""),dpi=300)

#ANOVA and Tukey's confirm the pulses have statistically equivalent x0 values

#group individuals
#cursed transformation to show all pulse lengths equidistant
# rheobase_sigmoidal_plot$pulse_length <- as.numeric(rheobase_sigmoidal_plot$pulse_length)
# 
# rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=x0, color=X),
#                                     show.legend=F) +
#   geom_point(show.legend = F) +
#   scale_x_continuous(labels = pulse_l, breaks = 1:7) +
#   geom_line(show.legend=F)
# 
# # Integrate MAMU into plots ----
# #commented out because i don't like these lol
# # rheobase_sigmoidal_plot %>% 
# #   ggplot(aes(MAMU, x0, color = as.factor(pulse_length)))+
# #   geom_point()
# # 
# # rheobase_sigmoidal_plot %>% 
# #   ggplot(aes(MAMU, x0))+
# #   geom_point() + 
# #   facet_wrap(~as.factor(pulse_length))
# 
# #pulse width vs EC10 -----
# rheobase_sigmoidal_plot$pulse_length <- as.factor(rheobase_sigmoidal_plot$pulse_length)
# 
# rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=EC10, fill=pulse_length)) +
#   geom_boxplot() +
#   scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#   geom_jitter(color="black", size=0.4, alpha=0.9) +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("Rheobase EC10 vs pulse width") +
#   xlab("Pulse width (us)")
# 
# #group individuals
# rheobase_sigmoidal_plot$pulse_length <- as.numeric(rheobase_sigmoidal_plot$pulse_length)
# 
# rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=EC10, color=X),
#                                     show.legend=F) +
#   geom_point(show.legend = F) +
#   scale_x_continuous(labels = pulse_l, breaks = 1:7) +
#   geom_line(show.legend=F)
# 
# #pulse width vs EC90 -----
# #some EC90 values are way too big, remove those
# rheobase_sigmoidal_plot$pulse_length <- as.factor(rheobase_sigmoidal_plot$pulse_length)
# 
# rheobase_sigmoidal_plot %>% 
#   filter(EC90<1000) %>% 
#   ggplot( aes(x=pulse_length, y=EC90, fill=pulse_length)) +
#   geom_boxplot() +
#   scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#   geom_jitter(color="black", size=0.4, alpha=0.9) +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("Rheobase EC90 vs pulse width") +
#   xlab("Pulse width (us)")
# 
# #group individuals
# rheobase_sigmoidal_plot$pulse_length <- as.numeric(rheobase_sigmoidal_plot$pulse_length)
# 
# rheobase_sigmoidal_plot %>% 
#   filter(EC90<1000) %>% 
#   ggplot( aes(x=pulse_length, y=EC90, color=X),
#                                     show.legend=F) +
#   geom_point(show.legend = F) +
#   scale_x_continuous(labels = pulse_l, breaks = 1:7) +
#   geom_line(show.legend=F)

#same as it previously was with spike at 10000 but no real trends

#IC50 plots (data formatted in analyses script for now) -----
# plot(x0 ~ IC50, data = IC50dat)
# 
# #plots for norm and bobby ----
# #storage dataframe
# all_rheo <- data.frame(
#   Dose = double(),
#   Ind = character(),
#   Snake = character(),
#   scaled_force = double(),
#   pulse_length = integer(),
#   MAMU = double()
# )
# 
# #another dataframe for the max force plot I'm going to make after this
# max_forces <- data.frame(
#   Ind = character(),
#   max_force = double(),
#   pulse_length = integer(),
#   MAMU = double(),
#   Snake = character()
# )
# 
# for (i in 1:length(pulse_l)) {
#   infile <- toString(paste("OutFiles/Rheobase/",foldername,"/Force_scaled/Rheo_",pulse_l[i],".csv", sep=""))
#   df <- read.csv(infile)
#   
#   #deal with max values first, to use later
#   max_long <- df[17,c(-1,-2)] %>% 
#     pivot_longer(
#       everything(),
#       names_to = "Ind",
#       values_to = "max_force"
#     ) %>% 
#     mutate(
#       Snake = sub("_.*", "", Ind),
#            pulse_length = pulse_l[i],
#            MAMU = 0
#     ) %>% 
#     rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") %>% #pull in MAMU
#   mutate(
#     pulse_length = as.factor(pulse_length)
#     )
#   max_forces <- rbind(max_forces,max_long)
#   
#   df <- df[-17,-1] #remove index line and maxes
#   dose_labels = df$Dose
#   df_long <- df %>%
#     pivot_longer(
#       cols = -Dose,
#       names_to = "Ind",
#       values_to = "scaled_force"
#     ) %>% 
#     mutate(
#       Snake = sub("_.*", "", Ind),
#       pulse_length = pulse_l[i],
#       MAMU = 0
#     ) %>% 
#     rows_update(distinct(keydf), by = "Snake", unmatched = "ignore") %>% #pull in MAMU
#     mutate(
#       pulse_length = as.factor(pulse_length),
#       Dose = as.numeric(Dose)
#     )
#   all_rheo <- rbind(all_rheo,df_long)
# }
# 
# for (page in 1:6) {
# plot <- ggplot(all_rheo, aes(x = Dose, y = scaled_force, group = pulse_length, color = pulse_length)) +
#   geom_point() + 
#   facet_wrap_paginate(~ Ind, ncol = 3, nrow = 2, page = page) +
#   geom_smooth(method = "nls",
#               method.args = list(formula = y ~  1 + (-1) / (1 + (x / x0)^dx),
#                                  start = list(x0 = 40 , dx = 1)
#               ),
#               data = all_rheo,
#               se = FALSE) +
#   scale_color_viridis(discrete = T) +
#   xlab("Current (mA)") + 
#   ylab("Proportion of maximal force")
# print(plot)
# }
# 
# #plot each set of raw data individually, one sheet per snake ----
# #i can't figure out how to make facet_wrap_paginate work well for odd numbers so i'm
# # adding a fake variable to each set just to give each snake 8 "pulse lengths"
# 
# #also they don't all have data for all plots
# 
# #best i can think of is removing all the 0,0 observations and then adding one
# # for each dataset (removing so there's no double data to mess up curve fitting)
# 
# no_0_rheo <- all_rheo[all_rheo$Dose != 0,]
# fake_base <- matrix(
#   c(
#     0,"Ind",0,"Snake",0,"MAMU",
#     0,"Ind",0,"Snake",50,"MAMU",
#     0,"Ind",0,"Snake",100,"MAMU",
#     0,"Ind",0,"Snake",200,"MAMU",
#     0,"Ind",0,"Snake",500,"MAMU",
#     0,"Ind",0,"Snake",1000,"MAMU",
#     0,"Ind",0,"Snake",10000,"MAMU",
#     0,"Ind",0,"Snake",50000,"MAMU"
#   ),
#   nrow= 8,
#   ncol = 6,
#   byrow=T,
#   
# )
# 
# colnames(fake_base) <- c("Dose","Ind","scaled_force","Snake","pulse_length","MAMU")
# fake_base<- as.data.frame(fake_base)
# 
# fake_zeros <- data.frame(
#   Dose = double(),
#   Ind = character(),
#   scaled_force = double(),
#   Snake = character(),
#   pulse_length = integer(),
#   MAMU = double()
# )
# 
# inds <- unique(all_rheo$Ind)
# 
# for(ind in inds){
#   fake_ind <- fake_base %>% 
#     mutate(
#       Dose = as.numeric(Dose),
#       Ind = ind,
#       scaled_force = as.numeric(scaled_force),
#       Snake =  sub("_.*", "", Ind),
#       MAMU = 0, #this should be fine for plotting purposes
#       pulse_length = as.factor(pulse_length)
#     )
#   fake_zeros <- rbind(fake_zeros,fake_ind)
# }
# 
# all_rheo_w_fake <- rbind(no_0_rheo,fake_zeros)
# all_rheo_w_fake <- all_rheo_w_fake[order(all_rheo_w_fake$Ind, all_rheo_w_fake$pulse_length,all_rheo_w_fake$Dose),]
# 
# #now make a page of plots for each snake
# for (page in 1:34) {
#   plot <- ggplot(all_rheo_w_fake, aes(x = Dose, y = scaled_force, 
#                                group = interaction(Ind,pulse_length), 
#                                color = pulse_length)) +
#     geom_point(
#       # aes(x = Dose, y = scaled_force, 
#       #             group = interaction(Snake,pulse_length))
#       ) + 
#     facet_wrap_paginate(~ Ind + pulse_length, ncol = 4, nrow = 2, page = page) +
#     geom_smooth(method = "nls",
#                 method.args = list(formula = y ~  1 + (-1) / (1 + (x / x0)^dx),
#                                    start = list(x0 = 60 , dx = 1)
#                                    # ,
#                                    # control = nls.lm.control(maxiter = 1024, maxfev = 1024),
#                                    # lower = c(0,-Inf), upper = c(600,Inf), algorithm = "port"
#                 ),
#                 data = all_rheo_w_fake, #don't want to fit one for the fake plot
#                 se = FALSE) +
#     scale_color_viridis(discrete = T) +
#     xlab("Current (mA)") + 
#     ylab("Proportion of maximal force")
#   print(plot)
#   ggsave(paste("OutFiles/Rheobase/",foldername,"/plots/rheo_",sort(unique(all_rheo$Ind))[page],".png",sep=""),dpi=300)
# }
# 
# #strength-duration curve with max force instead of 50%
# #general view
# 
# max_forces %>% ggplot( aes(x=pulse_length, y=max_force, fill=pulse_length)) +
#   geom_boxplot() +
#   scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#   geom_jitter(color="black", size=0.4, alpha=0.9) +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("Rheobase max force vs pulse width") +
#   xlab("Pulse width (us)")
# 
# #this is the reverse of what we expected as far as i can tell, it certainly
# # isn't the strength-duration curve we wanted
# 
# #group individuals
# max_forces %>% ggplot( aes(x=pulse_length, y=max_force, color=Ind, group = Ind),
#                                     show.legend=F) +
#   geom_point(show.legend = F) +
#   geom_line(show.legend=F)
# 
# # Integrate MAMU into plots ----
# max_forces %>% 
#   ggplot(aes(MAMU, max_force, color = as.factor(pulse_length)))+
#   geom_point()
# 
# #vaguely looks like contraction strength is weaker for the more resistant
# # snakes at higher pulse lengths, which i guess is neat?
# #need to run stats on this stuff now
# 
# max_forces %>% 
#   ggplot(aes(MAMU, max_force))+
#   geom_point() + 
#   facet_wrap(~as.factor(pulse_length))
# #okay not so much when it's pulled out like this

#loop through EC50 plots for each genotype
genotypes <- c("LVNV","WT_sirtalis","EPN","P","T","WT_elegans","WT_hammondii")
pulse_l <- c(50,100,200,500,1000,5000)

for (type in genotypes) {
  foldername = type
  
  rheobase_sigmoidal_plot <- read.csv(toString(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal_sub/Sigmoidal-rpt.csv",sep=""))) %>% 
    filter(pulse_length %in% c(50,100,200,500,1000,5000))
  
  #general view
  rheobase_sigmoidal_plot$pulse_length <- as.factor(rheobase_sigmoidal_plot$pulse_length)
  
  ec50 <- rheobase_sigmoidal_plot %>% ggplot( aes(x=pulse_length, y=x0, fill=pulse_length)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste(foldername,"Rheobase EC50 vs pulse width",sep=" ")) +
    xlab("Pulse width (us)")
  print(ec50)
  ggsave(paste("OutFiles/Rheobase/",foldername,"/plots/",foldername,"_EC50_sub200.png",sep=""),dpi=300)
  
  rheobase_sigmoidal_plot$pulse_length <- as.numeric(rheobase_sigmoidal_plot$pulse_length)
  
  ind_ec50 <- rheobase_sigmoidal_plot %>% 
    ggplot( aes(x=pulse_length, y=x0, color=X),
            show.legend=F) +
    geom_point(show.legend = F) +
    scale_x_continuous(labels = pulse_l, breaks = 1:6) +
    geom_line(show.legend=F) +
    ggtitle(paste(foldername,"Rheobase EC50 vs pulse width",sep=" "))
  print(ind_ec50)
  ggsave(paste("OutFiles/Rheobase/",foldername,"/plots/",foldername,"_EC50_inds_sub200.png",sep=""),dpi=300)
}

#average over genotypes -----
genotypes <- c("WT_sirtalis","LVNV","T")
geno_means <- matrix(ncol = 4)
colnames(geno_means) <- c("pulse_length","mean","se","genotype")

for (type in genotypes) {
  foldername = type
  rheo_ecs <- read.csv(paste("OutFiles/Rheobase/",foldername,"/Sigmoidal/Sigmoidal-rpt.csv",sep=""))
  
  rheo_ecs <- rheo_ecs[-1,-1] %>% 
    filter(MusMassg != 0.000999) %>% 
    group_by(pulse_length) %>% 
    summarize(
      mean = mean(x0), 
      se = sd(x0)/sqrt(length(x0))
    ) %>% 
    mutate(
      genotype = foldername
    )
  
  geno_means <- rbind(geno_means,rheo_ecs)
  
}
write.csv(geno_means[-1,],"OutFiles/Rheobase/All_geno_means.csv",row.names=F)

rheo_df <- read.csv("OutFiles/Rheobase/All_geno_means.csv")
rheo_df <- rheo_df %>% 
  filter(
    pulse_length %in% c(50,100,200,500,1000)
  ) %>% 
  mutate(
    genotype = as.factor(genotype)
  )
rheo_df$genotype <- factor(rheo_df$genotype, levels = c("WT_sirtalis","T","LVNV"))
rheo_df$pulse_length <- as.factor(rheo_df$pulse_length)

png(filename="OutFiles/Rheobase/plots/WT_rheoavg_dot.png",width=4500,height=2700,res=500)
WT_rheo <- ggplot(subset(rheo_df, genotype %in% c("WT_sirtalis","T")), 
                  aes(x=pulse_length, y=mean, colour=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 1,position = position_dodge(width=0.6),
                width = 0.5/3) +
  geom_line(linewidth = 1, linetype=2) +
  geom_point(size = 3,position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis")))),
                     values = c("#000000","green")) +
  scale_x_discrete(expand=c(0,1)) +
  scale_y_continuous(limits=c(0,150)) +
  #scale_x_continuous(breaks = c(50,100,200,500,1000)) +
  labs(y = "Current that Produces \n 50% Muscle Activation (mA)", 
       x = expression(paste("Pulse Width (",mu,"s)",sep="")),
       title = "Rheobase EC50 vs Pulse Width") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16),
        axis.text = element_text(size =14))

WT_rheo
dev.off()


png(filename="OutFiles/Rheobase/plots/WT_T_rheoavg_dot.png",width=4500,height=2700,res=500)
WT_T_rheo <- ggplot(subset(rheo_df, genotype %in% c("WT_sirtalis","T")), 
                    aes(x=pulse_length, y=mean, colour=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 1, position = position_dodge(width=0.6),
                width = 0.25) +
  geom_line(linewidth = 1,linetype=2) +
  geom_point(size = 3, position = position_dodge(width= 0.6)) +
  scale_x_discrete(expand=c(0,1)) +
  scale_y_continuous(limits=c(0,150)) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                                   expression(italic("Th. couchii"))),
                     values = c("#000000","#1C8642")) +
  #scale_x_continuous(breaks = c(50,100,200,500,1000)) +
  labs(y = "Current that Produces \n 50% Muscle Activation (mA)", 
       x = expression(paste("Pulse Width (",mu,"s)",sep="")),
       title = "Rheobase EC50 vs Pulse Width") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16),
        axis.text = element_text(size =14))

WT_T_rheo
dev.off()

png(filename="OutFiles/Rheobase/plots/WT_T_LVNV_rheoavg_dot.png",width=4500,height=2700,res=500)
WT_T_LVNV_rheo <- ggplot(subset(rheo_df, genotype %in% c("WT_sirtalis","T","LVNV")), 
                    aes(x=pulse_length, y=mean, colour=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 1, position = position_dodge(width=0.6),
                width = 0.5) +
  geom_line(linewidth = 1,linetype=2,aes(group = genotype)) +
  geom_point(size = 3,position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis"))),
                                                   expression(italic("Th. couchii")),
                                                   expression(paste("LVNV ",italic("Th. sirtalis")))),
                     values = c("#000000","#1C8642","#CC1E1D")) +
  scale_x_discrete(expand=c(0,1)) +
  scale_y_continuous(limits=c(0,150)) +
 # scale_x_continuous(breaks = c(50,100,200,500,1000)) +
  labs(y = "Current that Produces \n 50% Muscle Activation (mA)", 
       x = expression(paste("Pulse Width (",mu,"s)",sep="")),
       title = "Rheobase EC50 vs Pulse Width") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16),
        axis.text = element_text(size =14))

WT_T_LVNV_rheo
dev.off()

png(filename="OutFiles/Rheobase/plots/WT_rheoavg_line.png",width=4500,height=2700,res=500)
WT_rheo <- ggplot(subset(rheo_df, genotype %in% c("WT_sirtalis")), 
                  aes(x=pulse_length, y=mean, colour=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 1,position = position_dodge(width=0.6),
                width = 0.5/3) +
  geom_line(linewidth = 1, linetype=1,aes(group = genotype)) +
  geom_point(size = 3,position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("WT ",italic("Th. sirtalis")))),
                     values = c("#000000")) +
  scale_x_discrete(expand=c(0,1)) +
  scale_y_continuous(limits=c(0,150)) +
  #scale_x_continuous(breaks = c(50,100,200,500,1000)) +
  labs(y = "Current that Produces \n 50% Muscle Activation (mA)", 
       x = expression(paste("Pulse Width (",mu,"s)",sep="")),
       title = "Rheobase EC50 vs Pulse Width") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16),
        axis.text = element_text(size =14))

WT_rheo
dev.off()

png(filename="OutFiles/Rheobase/plots/T_rheoavg_line.png",width=4500,height=2700,res=500)
T_rheo <- ggplot(subset(rheo_df, genotype %in% c("T")), 
                  aes(x=pulse_length, y=mean, colour=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 1,position = position_dodge(width=0.6),
                width = 0.5/3) +
  geom_line(linewidth = 1, linetype=1,aes(group = genotype)) +
  geom_point(size = 3,position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste(italic("Th. couchii")))),
                     values = c("#1C8642")) +
  scale_x_discrete(expand=c(0,1)) +
  scale_y_continuous(limits=c(0,150)) +
  #scale_x_continuous(breaks = c(50,100,200,500,1000)) +
  labs(y = "Current that Produces \n 50% Muscle Activation (mA)", 
       x = expression(paste("Pulse Width (",mu,"s)",sep="")),
       title = "Rheobase EC50 vs Pulse Width") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16),
        axis.text = element_text(size =14))

T_rheo
dev.off()

png(filename="OutFiles/Rheobase/plots/LVNV_rheoavg_line.png",width=4500,height=2700,res=500)
LVNV_rheo <- ggplot(subset(rheo_df, genotype %in% c("LVNV")), 
                  aes(x=pulse_length, y=mean, colour=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 1,position = position_dodge(width=0.6),
                width = 0.5/3) +
  geom_line(linewidth = 1, linetype=1,aes(group = genotype)) +
  geom_point(size = 3,position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Genotype", labels = c(expression(paste("LVNV ",italic("Th. sirtalis")))),
                     values = c("#CC1E1D")) +
  scale_x_discrete(expand=c(0,1)) +
  scale_y_continuous(limits=c(0,150)) +
  #scale_x_continuous(breaks = c(50,100,200,500,1000)) +
  labs(y = "Current that Produces \n 50% Muscle Activation (mA)", 
       x = expression(paste("Pulse Width (",mu,"s)",sep="")),
       title = "Rheobase EC50 vs Pulse Width") +
  expand_limits(x = 0, y = 0) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 18),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=16),
        axis.text = element_text(size =14))

LVNV_rheo
dev.off()

#using the curve fitted with the rheobase equation -----
rheo <- read.csv("OutFiles/Rheobase/exp_decay_rheo_method.csv")

params_rheo <- rheo %>% group_by(Genotype) %>% 
  summarize(mean = mean(t),
            se = sd(t)/sqrt(length(t)))

x = seq(1,1000,1)
means <- matrix(nrow=length(x),ncol = nrow(params_rheo))


for (genotype in 1:nrow(params_rheo)) {
  rheobase = params_rheo$mean[genotype]
  rheo_se = params_rheo$se[genotype]
  y = rheobase*(1+(2*rheobase/x))
  
  means[,genotype] = y
}

colnames(means) <- params_rheo$Genotype
rownames(means) <- x

means <- as.data.frame(t(means)) 
 
means$Genotype = rownames(means)

plot_means <- pivot_longer(data = means, cols = -Genotype,
                 names_to = "pulse_width", values_to = "ec50") %>% 
  mutate(
    pulse_width = as.numeric(pulse_width)
  )

ggplot(plot_means %>% filter(Genotype != "P")) +
  geom_line(aes(x = pulse_width,y=ec50,group=Genotype, color=Genotype),
            lwd = 1) + xlim(0,200) + ylim(0,500)

ggsave("Plots/rheobase_plots/strength_duration.png",dpi=300)
#+  ylim(0,170)

#now do it again with exponential decay function ----
#this equation does not fit well
# tau <- read.csv("OutFiles/Rheobase/exp_decay_method2.csv")
# 
# params_exp <- tau %>% group_by(Genotype) %>% 
#   summarize(mean_t = mean(t),
#             se_t = sd(t)/sqrt(length(t)),
#             mean_x0 = mean(x0),
#             se_x0 = sd(x0)/sqrt(length(x0)))
# 
# x = seq(1,10000,1)
# means_exp <- matrix(nrow=length(x),ncol = nrow(params_exp))
# 
# for (genotype in 1:nrow(params_exp)) {
#   x0 = params_exp$mean_x0[genotype]
#   t = params_exp$mean_t[genotype]
#   
#   y = x0*exp(-x*t)
#   means_exp[,genotype] = y
# }
# 
# colnames(means_exp) <- params_exp$Genotype
# rownames(means_exp) <- x
# 
# means_exp <- as.data.frame(t(means_exp)) 
# 
# means_exp$Genotype = rownames(means_exp)
# 
# plot_means <- pivot_longer(data = means_exp, cols = -Genotype,
#                            names_to = "pulse_width", values_to = "ec50") %>% 
#   mutate(
#     pulse_width = as.numeric(pulse_width)
#   )
# 
# ggplot(plot_means) +
#   geom_line(aes(x = pulse_width,y=ec50,group=Genotype, color=Genotype)) + xlim(0,500) + ylim(0,200)
