#This script consolidates, manipulates, and explores the organismal-level
# data (primarily sprint speed and morphological data); most steps were used
# to prepare the data for the models run in script 2 and to determine which models 
# would be the most useful/may show a relationship

#need to move the working directory because the default for this R project is 
# in the Myography folder but I don't want to make a new project just for this
setwd('../tradeoff_modeling')
getwd()

#Reading in data ----
couchiidat <- read.csv("data_raw/couchiionlydata.csv")
newdat <- read.csv("data_raw/summer24couchii.csv")

# Libraries: -------
library(tidyverse)
library(inspectdf)
library(ggplot2)
library(DHARMa)
library(SuppDists)
library(lme4)
library(sjPlot)
library(lavaan)
library(lavaanPlot)

#1: Data formatting ----
couchiidat <- couchiidat %>% 
  bind_rows(newdat) %>% 
  select(species, INDIV, CollPreNo, locality, County, mass_g, Sex,
         SVL_cm, TL_cm, baseline_ms, MAMU) %>% 
  mutate ( #misc small changes
    MAMU = as.double(MAMU),
    species = as.factor(species),
    log_MAMU = log(MAMU),
    Sex = as.factor(Sex)
  )

#Making a scaled data-frame for modeling:
scaledat <- couchiidat %>% 
  select(species, INDIV, mass_g, SVL_cm, TL_cm, baseline_ms, MAMU, log_MAMU, Sex) %>% 
  mutate(
    mass_g = scale(mass_g, center = F)[,1],
    SVL_cm = scale(SVL_cm, center = F)[,1],
    TL_cm = scale(TL_cm, center = F)[,1],
    MAMU = scale(MAMU, center = F)[,1],
    log_MAMU = scale(log_MAMU, center = F)[,1] 
    #feels odd to scale something that's already logged, but I'll compare 
    #using it to using regular MAMU scaled (because it would also feel wrong
    #to scale all but one variable)
  )

#2: Data exploration -----
#2A: Looking at data itself -----

#examine dataset:
summary(couchiidat)
inspect_types(couchiidat)

datinfo <- inspect_cat(couchiidat)
show_plot(datinfo)
#I don't expect to be able to use any of these for modeling, maybe locality once
# I go back through it later

par(mfrow = c(2,1))
hist(couchiidat$MAMU)
#this is definitely going to cause some problems, I'll experiment with log-
# transforming it in some of the models
hist(couchiidat$log_MAMU)
#that's much more workable, even if it'll make interpretation slightly harder
hist(couchiidat$mass_g)
#squite skewed
hist(couchiidat$SVL_cm)
# not quite normal but not too skewed
hist(couchiidat$baseline_ms)
#also about as expected

#2B: Plotting for possible interactions ----
# first, look at relationship between MAMU and speed 
ggplot(couchiidat, aes(MAMU, baseline_ms))+
  geom_point() + geom_smooth(method='lm')
ggplot(couchiidat, aes(log_MAMU, baseline_ms))+
  geom_point() + geom_smooth(method='lm')
#both of these show lovely flat lines

ggplot(couchiidat, aes(SVL_cm, baseline_ms))+
  geom_point() + geom_smooth(method='lm')

#I know Chris wants to keep those fast neonates but i don't know how we'll justify
# the negative association of SVL and speed - maybe motivation?

ggplot(couchiidat, aes(mass_g, baseline_ms))+
  geom_point() + geom_smooth(method='lm')

#same thing, weird relationship here