
#Reading in data ----
couchiidat <- read.csv("tradeoff_modeling/couchiionlydata.csv")
newdat <- read.csv("tradeoff_modeling/summer24couchii.csv")

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

#3: GLMs and GLMMs -----
#3A: Selecting distribution and link function using GLMs -----
# I don't know if this is theoretically supported, but I wanted to start by
# selecting distribution and link functions with simpler GLMs before adding in
# random effects

#3Ai: Normal distribution (for baseline): ----
glm_norm <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                family = gaussian(link = "identity"), data = scaledat)
summary(glm_norm)
#some diagnostics with DHARMa (because it'll be good for the non-normal dists)
simulateResiduals(glm_norm, plot = T)
#well, the QQ residuals actually look good but not the quantiles of the residuals

plot_model(glm_norm, show.values = T, show.p = F)

#Normal with log link (for curiosity):
glm_norm_log <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                    family = gaussian(link = "log"), data = scaledat)
summary(glm_norm_log)
simulateResiduals(glm_norm_log, plot = T)
#quantile plot is worse, residual deviance slightly closer to null deviance and 
# therefore worse but not terribly so

#3Aii: Gamma distribution: -----
#First I'll check the default gamma link function, the inverse
glm_gamma_inv <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                     family = Gamma(link = "inverse"), data = scaledat)
summary(glm_gamma_inv)
#based on AIC alone, this model is better than the normal model
#null-residual comparison is about the same though
#however it does say mamu has a positive effect on speed which we know shouldn't be true

glm_gamma_inv_res <- simulateResiduals(glm_gamma_inv, plot = T)
testResiduals(glm_gamma_inv_res)
#uniformity looks pretty good, residuals aren't over- or under-dispersed, only
# a small selection of outliers; still fails quantile test
inv <- function(x) {
  1/x
}

plot_model(glm_gamma_inv, type = "std", sort.est = T, show.values = T,
           title = "Coefficients from model with inverse", 
           transform = "inv")
#I'm not sure i actually got this effect estimating working correctly, I need to 
# check that somehow 

#now the log link function, which I read is more commonly used
glm_gamma_log <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                     family = Gamma(link = "log"), data = scaledat)
summary(glm_gamma_log)
#according to AIC, this is an even better model

glm_gamma_log_res <- simulateResiduals(glm_gamma_log, plot = T)
testResiduals(glm_gamma_log_res)
#again, quantile deviation but otherwise okay

plot_model(glm_gamma_log, type = "est", show.values = T,
           title = "Coefficients from model with log link", show.p = F,
           transform = "exp")
#more sensible results here

#3Aiii: Inverse Gaussian: ----
#The gamma distribution wasn't a bad fit, but it could clearly be improved upon,
# so we'll check out the inverse gaussian distribution last
#first just using the default link function, 1/mu^2
glm_invg_mu <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                   family = inverse.gaussian(link = "1/mu^2"), data = scaledat)
summary(glm_invg_mu)
#Based on both AIC and deviance, this model is worse than both gamma fits
#I'll check using the residuals again:
glm_invg_mu_res <- simulateResiduals(glm_invg_mu, plot = T)
testResiduals(glm_invg_mu_res)
#quantiles fail yet again

#maybe a different link function will help:
glm_invg_inv <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                    family = inverse.gaussian(link = "inverse"), data = scaledat)
summary(glm_invg_inv)
#both inverse gaussians worse than either gamma distribution model
glm_invg_inv_res <- simulateResiduals(glm_invg_inv, plot = T)
testResiduals(glm_invg_inv_res)
#same quantile issues again

#one last try, the log link function:
glm_invg_log <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g, 
                    family = inverse.gaussian(link = "log"), data = scaledat)
summary(glm_invg_log)
#Based on AIC alone this model is not better
glm_invg_log_res <- simulateResiduals(glm_invg_log, plot = T)
testResiduals(glm_invg_log_res)
#more errors this time, and failed quantile thing again

#gamma with log link seems the best for this dataset, based just on AIC since 
# they all have quantile issues

#3B: GLMMs ----
#I nuked these because the dataset is currently empty of things that would make
# sense as random effects; once i fix the locality section i'll do these

#4: Path models -----
path1 = '
  baseline_ms ~ log_MAMU 
  baseline_ms ~ SVL_cm 
  baseline_ms ~ mass_g'

path1.fit = sem(path1, data = scaledat)
lavResiduals(path1.fit)
summary(path1.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "path1", model = path1.fit, coefs = TRUE)
#not great fit based on diagnostics

path2 = '
  mass_g ~ SVL_cm
  baseline_ms ~ log_MAMU + SVL_cm + mass_g'

path2.fit <- sem(path2, data = scaledat)
lavResiduals(path2.fit)
summary(path2.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "path2", model = path2.fit, coefs = TRUE)
#also showing that negative effect of svl that i just don't believe
#it is neat to confirm that the SVL~mass effect doesn't affect the other paths

#tbh i'm a little stumped on what else to do besides fixing the localities for
# a glmm and the quantile thing (and those weird data points), so much of the
# fun in my research design analyses came from the multispecies aspect