#This script makes GLMs (and eventually GLMMs) and path models to explore 
# the effects of various variables, including MAMU, on sprint speed

#need to move the working directory because the default for this R project is 
# in the Myography folder but I don't want to make a new project just for this
setwd('../tradeoff_modeling')
getwd()

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

#GLMs and GLMMs -----
#1: Selecting distribution and link function using GLMs -----
# I don't know if this is theoretically supported, but I wanted to start by
# selecting distribution and link functions with simpler GLMs before adding in
# random effects

#1A: Normal distribution (for baseline): ----
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

#1B: Gamma distribution: -----
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

#1C: Inverse Gaussian: ----
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

#2: GLMMs ----
#I nuked these because the dataset is currently empty of things that would make
# sense as random effects; once i fix the locality section i'll do these

#3: Path models -----
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