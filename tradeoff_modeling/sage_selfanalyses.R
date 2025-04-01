# Analyses of multispeciesdata.csv for Research Design ----
# Includes: data exploration, data cleaning/manipulation (originally from script
# sage_datacleaning.R), mixed models to explore distributions, basic SEM to
# (attempt to) get better idea of what relationships to include, and modeling
# of MAMU from alleles and mutations to compare to literature results

#multidat <- read.csv("multispeciesdata.csv")

# Libraries: -------
library(tidyverse)
library(inspectdf)
library(ggplot2)
library(DHARMa)
library(SuppDists)
library(lme4)
library(lavaan)
library(lavaanPlot)

#1: Data formatting ------
#I'm making columns for each allele and its count in each individual
#0 = not there, 1 = heterozygous, 2 = homozygous
#V2 coded directly from LVNV because it's only found in that combination,
#T also coded directly because of issues with LVNV, otherwise all found with grepl

#any suggestions on how to make the dummy variable step cleaner would be 
# appreciated, I looked up some shorter tidyverse syntax but kept getting errors
# I couldn't figure out how to reconcile
multidat <- multidat %>% 
  mutate ( #misc small changes
    MAMU = as.double(MAMU),
    species = as.factor(species),
    region = as.factor(region),
    age = as.factor(age),
    log_MAMU = log(MAMU)
  ) %>% 
  mutate( #allele a
    E_allele = ifelse(grepl("E",allele_a,fixed=T),1,0),
    P_allele = ifelse(grepl("P",allele_a,fixed=T),1,0),
    N_allele = ifelse(grepl("N",allele_a,fixed=T),1,0),
    T_allele = ifelse(allele_a == "T",1,0),
    L_allele = ifelse(grepl("L",allele_a,fixed=T),1,0),
    V1_allele = ifelse(grepl("V",allele_a,fixed=T),1,0),
    V2_allele = ifelse(grepl("LVNV",allele_a,fixed=T),1,0),
    A_allele = ifelse(grepl("A",allele_a,fixed=T),1,0),
  ) %>% 
  mutate ( #add on allele b
    E_allele = ifelse(grepl("E",allele_b,fixed=T),E_allele+1,E_allele),
    P_allele = ifelse(grepl("P",allele_b,fixed=T),P_allele+1,P_allele),
    N_allele = ifelse(grepl("N",allele_b,fixed=T),N_allele+1,N_allele),
    T_allele = ifelse(allele_b == "T",T_allele+1,T_allele),
    L_allele = ifelse(grepl("L",allele_b,fixed=T),L_allele+1,L_allele),
    V1_allele = ifelse(grepl("V",allele_b,fixed=T),V1_allele+1,V1_allele),
    V2_allele = ifelse(grepl("LVNV",allele_b,fixed=T),V2_allele+1,V2_allele),
    A_allele = ifelse(grepl("A",allele_b,fixed=T),A_allele+1,A_allele),
  ) 

#Making an alleles and MAMU (and species)-only data frame for part 5:
alleledat <- multidat %>% 
 # select(INDIV, species, MAMU, log_MAMU, allele_a, allele_b, E_allele:A_allele) %>% 
  subset(!is.na(MAMU)) %>% 
  subset(!is.na(allele_a)) %>% 
  subset(!is.na(allele_b)) %>% 
  mutate(ones = 1, fill = 0) %>% 
  pivot_wider(names_from = allele_b, values_from = ones, values_fill = 0) %>% 
  mutate(
    WT = ifelse(allele_a == "WT", WT + 1, WT),
    P = ifelse(allele_a == "P", P + 1,P),
    EPN = ifelse(allele_a == "EPN", EPN + 1, EPN),
    T = ifelse(allele_a == "T", T + 1, T),
    LVNV = ifelse(allele_a == "LVNV", LVNV + 1, LVNV),
    VA = ifelse(allele_a == "VA", VA + 1, VA),
    V = ifelse(allele_a == "V", V + 1, V)
  ) 


#Making a scaled data-frame for modeling:
scaledat <- multidat %>% 
  select(species, INDIV, region, age, mass_g, SVL_cm, baseline_ms, MAMU, log_MAMU,
         E_allele:A_allele) %>% 
  mutate(
    mass_g = scale(mass_g, center = F)[,1],
    SVL_cm = scale(SVL_cm, center = F)[,1],
    MAMU = scale(MAMU, center = F)[,1],
    log_MAMU = scale(log_MAMU, center = F)[,1] 
    #feels odd to scale something that's already logged, but I'll compare 
    #using it to using regular MAMU scaled (because it would also feel wrong
    #to scale all but one variable)
  )
# I didn't originally make this, but the GLMMs were upset about the scale of the
# data (or lack of scaling thereof) so here's the answer to that

#2: Data exploration -----
#2A: Looking at data itself -----
#(this and section 1 fed into each other, they just got reorganized to sort out
# the different parts of the workflow)

#examine dataset:
summary(multidat)
inspect_types(multidat)

datinfo <- inspect_cat(multidat)
show_plot(datinfo)
#age, alleles, locality, and species have the most completed data
# but locality is so variable that region is more likely to be useful
#sex basically empty and probably not imputable from the lack of data
#age possibly imputable because it's so strongly size-correlated

par(mfrow = c(2,1))
hist(multidat$MAMU)
#this is definitely going to cause some problems, I'll experiment with log-
# transforming it in some of the models
hist(multidat$log_MAMU)
#that's much more workable, even if it'll make interpretation slightly harder
hist(multidat$mass_g)
#about as expected
hist(multidat$SVL_cm)
# neat to see the spike of juveniles and otherwise semi-normal distribution
hist(multidat$baseline_ms)
#also about as expected

#2B: Plotting for possible interactions ----
# first, look at relationship between MAMU and speed by and w/o species
ggplot(multidat, aes(log_MAMU, baseline_ms))+
  geom_point() + geom_smooth()
ggplot(multidat, aes(log_MAMU, baseline_ms, color = species))+
  geom_point() + geom_smooth()
#these are both weird and the clean relationship predicted by theory isn't there
#but, it does appear that there's something species-specific here

ggplot(multidat, aes(SVL_cm, baseline_ms))+
  geom_point() + geom_smooth()
ggplot(multidat, aes(SVL_cm, baseline_ms, color = species))+
  geom_point() + geom_smooth()
#I don't know what's going on with couchii here, it almost feels like speed was
# measured in different units but it'd be unrealistically small otherwise

ggplot(multidat, aes(mass_g, baseline_ms))+
  geom_point() + geom_smooth()
ggplot(multidat, aes(mass_g, baseline_ms, color = species))+
  geom_point() + geom_smooth()
#that chunk of fast couchii neonates really throws off what would otherwise be
# a reasonable logarithmic curve
ggplot(multidat, aes(mass_g, baseline_ms, color = age))+
  geom_point() + geom_smooth()
# and this makes me inclined to think that something is wrong because none of 
# the other neonates are nearly that fast

#I wonder if they might have been recorded in feet/s for that group (the data
# sheet didn't have units in it) but there's no metadata so I can't prove anything

#I'm going to just leave it for now and revisit it in the future, maybe see if
# I can figure out if anyone remembers how they recorded data in 2001 (doubtful)

#3: GLMs and GLMMs -----
#I'm going to try making GLMMs with both the Gamma and inverse gaussian
# distributions, and test a variety of link functions (because the log link
# is apparently most common, but I'm not sure the multiplicative effects it
# assumes the predictors have are appropriate for the system)

#3A: Selecting distribution and link function using GLMs -----
# I don't know if this is theoretically supported, but I wanted to start by
# selecting distribution and link functions with simpler GLMs before adding in
# random effects

#3Ai: Normal distribution (for baseline): ----
glm_norm <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                family = gaussian(link = "identity"), data = scaledat)
summary(glm_norm)
#some diagnostics with DHARMa (because it'll be good for the non-normal dists)
simulateResiduals(glm_norm, plot = T)
#well, the QQ residuals actually look good but the quantiles of the residuals
# definitely do not
plot_model(glm_norm, show.values = T, show.p = F)

#Normal with log link (for curiosity):
glm_norm_log <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                family = gaussian(link = "log"), data = scaledat)
summary(glm_norm_log)
simulateResiduals(glm_norm_log, plot = T)
#this is seemingly worse, KS test now says this doesn't come from a normal dist
# (which we know because the normal has negative support but still)

#3Aii: Gamma distribution: -----
#First I'll check the default gamma link function, the inverse
glm_gamma_inv <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                    family = Gamma(link = "inverse"), data = scaledat)
summary(glm_gamma_inv)
#based on AIC alone, this model is better than the normal model
#and null deviance higher than residual deviance, implying improved fit 
# when compared to intercept-only model
glm_gamma_inv_res <- simulateResiduals(glm_gamma_inv, plot = T)
testResiduals(glm_gamma_inv_res)
#uniformity looks pretty good, residuals aren't over- or under-dispersed, only
# a small selection of outliers
plot_model(glm_gamma_inv, type = "std", sort.est = T, show.values = T,
           title = "Coefficients from model with inverse", 
           transform = "inv")

inv <- function(x) {
  1/x
}
#now the log link function, which I read is more commonly used
glm_gamma_log <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                     family = Gamma(link = "log"), data = scaledat)
summary(glm_gamma_log)
#according to AIC, this is an even better model
#also, residual deviance slightly dec but null deviance about the same
glm_gamma_log_res <- simulateResiduals(glm_gamma_log, plot = T)
testResiduals(glm_gamma_log_res)
#though the model fits better, there is now deviance from the expected dispersion
# of residuals between fitted and simulated, and observations appear over-
# dispersed compared to the model we're using
plot_model(glm_gamma_log, type = "est", show.values = T,
           title = "Coefficients from model with log link", show.p = F,
           transform = "exp")

#3Aiii: Inverse Gaussian: ----
#The gamma distribution wasn't a bad fit, but it could clearly be improved upon,
# so we'll check out the inverse gaussian distribution last
#first just using the default link function, 1/mu^2
glm_invg_mu <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                     family = inverse.gaussian(link = "1/mu^2"), data = scaledat)
summary(glm_invg_mu)
#Based on both AIC and deviance, this model is worse than both gamma fits
#I'll check using the residuals again:
glm_invg_mu_res <- simulateResiduals(glm_invg_mu, plot = T)
testResiduals(glm_invg_mu_res)
#no overdispersion issue this time, but a clear issue with the fit to the dist

#maybe a different link function will help:
glm_invg_inv <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                   family = inverse.gaussian(link = "inverse"), data = scaledat)
summary(glm_invg_inv)
#Better than the last inverse gaussian based on AIC, but the deviance is still
# worse than either gamma distribution model
glm_invg_inv_res <- simulateResiduals(glm_invg_inv, plot = T)
testResiduals(glm_invg_inv_res)
#same issue with the KS test results as last time

#one last try, the log link function:
glm_invg_log <- glm(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species, 
                   family = inverse.gaussian(link = "log"), data = scaledat)
summary(glm_invg_log)
#Based on AIC alone this model is way better, and the difference between null
# and residual deviance is pretty good
glm_invg_log_res <- simulateResiduals(glm_invg_log, plot = T)
testResiduals(glm_invg_log_res)
#still, big issues with the dispersion test (though visually the QQ plot looks 
# good, I have to say)

#I think for the mixed models I'll use the gamma distribution with the inverse
# link function (only one that didn't have issues with KS or dispersion) and
# the inverse gaussian with a log link (best AIC despite KS test results showing
# it isn't actually the generating distribution)

#3B: GLMMs ----
#3Bi: Gamma distribution-based models 
#first I'll just add a random intercept by site
gamma_glmm1 <- glmer(baseline_ms ~ MAMU + SVL_cm*mass_g + species + (1|region),
                     data = scaledat, family = Gamma(link = "inverse"))
summary(gamma_glmm1)
testResiduals(simulateResiduals(gamma_glmm1))
#don't understand how the Gamma was good based on residuals before adding a
# random effect but was horribly worse afterwards

#try adding species as a random effect (justifiable depending on how my question
# is framed, and I'm curious what it'll do to the residuals)
gamma_glmm2 <- glmer(baseline_ms ~ MAMU + SVL_cm*mass_g + (1|species) ,
                     data = scaledat, family = Gamma(link = "inverse"))
summary(gamma_glmm2)
#this one didn't converge, so definitely not any better

#3Bii: inverse gaussian distribution-based models
#Trying the other distribution that worked semi-well:
invg_glmm1 <- glmer(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species + (1|region), 
                    family = inverse.gaussian(link = "log"), data = scaledat)
summary(invg_glmm1)
#this one also didn't converge so it's even worse; no inverse gaussian here

#3Biii: back to a normal distribution (just to see what happens)
norm_glmm <- glmer(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species + (1|region), 
                   family = gaussian(link = "identity"), data = scaledat)
summary(norm_glmm)
testResiduals(simulateResiduals(norm_glmm))
#there's negative skew to the data according to the QQ plot
#trying a log link function:
norm_glmm2 <- glmer(baseline_ms ~ log_MAMU + SVL_cm*mass_g + species + (1|region), 
                   family = gaussian(link = "log"), data = scaledat)
summary(norm_glmm2)
testResiduals(simulateResiduals(norm_glmm2))
#doesn't seem like changing the link function fixes the skew

# One last mixed model try: overhauling the model statement entirely to the one
# suggested in the suggested analyses document
gamma_glmm3 <- glmer(baseline_ms ~ log_MAMU + E_allele + P_allele + N_allele +
                      T_allele + L_allele + V1_allele + V2_allele + A_allele +
                       (1 | region),
                     family = Gamma(link = "inverse"), data = scaledat)
summary(gamma_glmm3)
#didn't converge and had to drop fixed effects because it was rank deficient

#TLDR: I would really appreciate feedback on what might have been happening with
# the fit of the gamma distribution changing drastically as soon as a random 
# effect was added, I've never seen anything like that and can't find anything
# about it online. For the presentation I'm going to focus on the GLM results
# and see if anyone in my group had more success with models than me, but long-
# term I want to know what was happening here

#4: SEM ------
#After that mess with GLMMs, let's see if the more complex SEMs we can make
# can capture the patterns we observe in our data (and compare different
# hypotheses about relationships between variables)

#4A: Simplest model ----
# The simplest model has SVL, mass, and MAMU as predictors of speed, with SVL
# and mass interacting to produce an endogenous body condition variable that
# is also predictive

model1 <- '
  # regression
  baseline_ms ~ bodcon + SVL_cm + mass_g + MAMU
  
  # latent variable definition
  bodcon =~ SVL_cm + mass_g
'
#run model
model1.fit <- sem(model1, data = multidat)
#the error messages make it clear this model isn't correct (and therefore can't
# be fit, giving the weird negative variances)

lavaanPlot(model1.fit, coefs = T)

#return Rsq value, standardize data, and estimate model fit (for later comparison)
summary(model1.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE) 

#I want to try the same model but with the MAMU logged:
model1.5 <- '
  # regression
  baseline_ms ~ bodcon + SVL_cm + mass_g + log_MAMU
  
  # latent variable definition
  bodcon =~ SVL_cm + mass_g
'
#run model
model15.fit <- sem(model1.5, data = multidat)
#well, that didn't help very much

#Maybe having SVL and mass predict both body condition and speed is confounding,
# so I'll try the actual simplest model (that still merits SEM):

model0 <- '
  # regression
  baseline_ms ~ bodcon + log_MAMU
  
  # latent variable definition
  bodcon =~ SVL_cm + mass_g
'

model0.fit <- sem(model0, data = multidat)
#this one at least tried to run and gave one less warning
summary(model0.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE) 
#still doesn't look like a very good model

#I noticed from the plots that a lot of the couchii observations don't follow
# the trend set by the other three species, so what if we remove them and try the
# model again? 

ncdat <- multidat[multidat$species != "Thamnophis couchii",]

model1.fit2 <- sem(model1, data = ncdat)
#no, this didn't help with all the errors either
model0.fit2 <- sem(model0, data = ncdat)
#this was also marginally better than the original try but not truly functional

#4B: Trying a more complex model ----
# I want to see if a more complex generative model helps at all (unlikely but
# I am really trying here)

#This one will dig into the TTX resistance side a bit more
modelttx <- '
  baseline_ms ~ SVL_cm + mass_g + TTXres
  TTXres ~ E_allele + P_allele + N_allele + L_allele + V1_allele + 
            V2_allele + A_allele
  TTXres =~ log_MAMU
'

#excluding couchii data from start because it lacks region data and all has the
# same alleles

modelttx.fit <- sem(modelttx, data = ncdat)
#no go on this either, time to just try things not in my suggested analyses
# because I don't think this relationship really can be modeled (which tracks
# with some of the papers I've read, half estimate an effect and half don't)

#5: MAMU vs alleles present -----
#As an alternative, I want to try to model MAMU as a result of the alleles
# present and then compare the estimated effect sizes from the model to the
# effect sizes from the literature on these alleles
#I'll do this both including and excluding couchii, because I think they'll 
# largely be confounding given the lack of genetic variation but I'm curious
# what effect they'll have on the estimated means and variance

#5A: Including Th. couchii ----
#I also want to make two models of each type, one with individual mutation 
# counts and one with mutation counts (the first is more interpretable but the
# second is more easily compared to the literature)
#allele glm:
allele_glm <- glm(MAMU ~ N + P + EPN + WT + EP + T + LVNV + VA + V, 
               data = alleledat)
summary(allele_glm)
#says it can't estimate an intercept for V, supposedly because of correlation
# between values (which I guess makes sense, most will be 0 simultaneously)
#unsurprisingly, EPN and LVNV have large estimated effect sizes

#test residuals to see if gamma dist would be better:
testResiduals(simulateResiduals(allele_glm))
#big diversion from the residuals, let's try a gamma
allele_glm2 <- glm(MAMU ~ N + P + EPN + WT + EP + T + LVNV + VA + V, 
                  data = alleledat, family = Gamma(link = "inverse"))
summary(allele_glm2)
#still couldn't est for V, has more significant coef ests but worse AIC
testResiduals(simulateResiduals(allele_glm2))
#despite AIC, residuals look way better -but- dispersion test is worse

#test another link function:
allele_glm3 <- glm(MAMU ~ N + P + EPN + WT + EP + T + LVNV + VA + V, 
                   data = alleledat, family = Gamma(link = "log"))
summary(allele_glm3)
#still couldn't est for V
testResiduals(simulateResiduals(allele_glm3))
#residuals look less bad according to KS test but only slightly; dispersion test 
# still bad, don't have to worry about whether multiplicative effects are 
# supported by theory for all the homozygotes but possibly wrong for heterozygotes

#I'll compare results of all three models to the literature-based effect sizes
# in my presentation

#individual mutation glm:
mut_glm <- glm(MAMU ~ E_allele + P_allele + N_allele + T_allele + L_allele + 
                 V1_allele + V2_allele + A_allele, data = alleledat )
summary(mut_glm)
#this one wasn't able to estimate the effect of V2, but the covariance makes
# more sense here because it's always found with L, V1, and N
#deviance is insane on this model though
testResiduals(simulateResiduals(mut_glm))
#same shaped QQ plot as last one, I presume the gamma will help again

#gamma with inverse link (default):
mut_glm2 <- glm(MAMU ~ E_allele + P_allele + N_allele + T_allele + L_allele + 
                 V1_allele + V2_allele + A_allele, data = alleledat,
                family = Gamma(link = "inverse"))
summary(mut_glm2)
#same issue with V2, but much better deviance (and AIC)
testResiduals(simulateResiduals(mut_glm2))
#and same issues with gamma and inverse link as last time, minor dispersion on
# QQ plot and overdispersion (though overdispersion isn't surprising when we 
# know the alleles don't account for all variation in MAMU)

#I'll also do a log model because some literature suggests that the mutations
# are (somehow) multiplicative in effect (and some lit suggests they're additive)
mut_glm3 <- glm(MAMU ~ E_allele + P_allele + N_allele + T_allele + L_allele + 
                  V1_allele + V2_allele + A_allele, data = alleledat,
                family = Gamma(link = "log"))
summary(mut_glm3)
#same issue with V2; residual deviance slightly worse which suggests mutations
# are more additive than multiplicative in effect
testResiduals(simulateResiduals(mut_glm3))
#still overdispersed but less so than the last one, still, the residuals imply
# the inverse link is probably better


#5B: Excluding Th. couchii -----
#need version of allele data frame without couchii
ncdat2 <- alleledat[alleledat$species != "Thamnophis couchii",]

#run same models as above (going straight to gamma distribution)
#allele-based model:
allele_glm_nc <- glm(MAMU ~ N + P + EPN + WT + EP + LVNV + VA + V, 
                   data = alleledat, family = Gamma(link = "inverse"))
summary(allele_glm_nc)
#able to estimate for each term! finally
testResiduals(simulateResiduals(allele_glm_nc))
#issues with KS test and overdispersion, but overdispersion honestly isn't as
# bad as a lot of the other models

#check out the log link again:
allele_glm_nc2 <- glm(MAMU ~ N + P + EPN + WT + EP + LVNV + VA + V, 
                     data = alleledat, family = Gamma(link = "log"))
summary(allele_glm_nc2)
#able to estimate for each term again; residual deviance slightly better but
# AIC worse
testResiduals(simulateResiduals(allele_glm_nc2))
#KS test better but overdispersion worse

#mutation-based model:
mut_glm_nc <- glm(MAMU ~ E_allele + P_allele + N_allele + L_allele + 
                  V1_allele + V2_allele + A_allele, data = ncdat2,
                family = Gamma(link = "inverse"))
summary(mut_glm_nc)
#still can't estimate for the V2 allele, but residual deviance and AIC both went
# way down which is good (and logical)
testResiduals(simulateResiduals(mut_glm_nc))
#however, more dispersion and an issue with the KS test
#data seems slightly right-skewed, which makes some sense because there's now
# a lower resistant:non-resistant snake ratio and therefore overall lower MAMUs

#try a log link in case taking out couchii changed anything:
mut_glm_nc2 <- glm(MAMU ~ E_allele + P_allele + N_allele + L_allele + 
                    V1_allele + V2_allele + A_allele, data = ncdat2,
                  family = Gamma(link = "log"))
summary(mut_glm_nc2)
#still can't estimate for the V2 allele, but residual deviance and AIC both went
# way down which is good (and logical)
testResiduals(simulateResiduals(mut_glm_nc2))
#diagnostics basically identical, residuals messed up in a slightly different 
# spot and the overdispersion is worse even though residual deviance is lower

#don't think we can conclude from either of these models whether mutation/allele
# effects are more additive or multiplicative, I think it may vary between cases
# (and I think if we had all heterozygotes we might be able to estimate better
# for the allele case at least)

#5C: Comparing above models with McFadden's R-squared ----
#I read that this metric can be used to estimate fit of glms in R; please let 
# me know if it isn't actually used/has fallen out of favor/is biased in some
# way

#not calculating for the normal dist-based models because of their poor fits

#allele-based:
#+couchii, inverse link
with(summary(allele_glm2), 1 - deviance/null.deviance)
#0.5463518
#+couchii, log link
with(summary(allele_glm3), 1 - deviance/null.deviance)
#0.5709613
#-couchii, inverse link
with(summary(allele_glm_nc), 1 - deviance/null.deviance)
#0.5465449
#-couchii, log link
with(summary(allele_glm_nc2), 1 - deviance/null.deviance)
#0.5709613

#mutation based:
#+couchii, inverse link
with(summary(mut_glm2), 1 - deviance/null.deviance)
#0.5463518
#+couchii, log link
with(summary(mut_glm3), 1 - deviance/null.deviance)
#0.5706402
#-couchii, inverse link
with(summary(mut_glm_nc), 1 - deviance/null.deviance)
#0.6654379
#-couchii, log link
with(summary(mut_glm_nc2), 1 - deviance/null.deviance)
#0.6967245 

#all very good R-squared values all things considered, and right about the 
# percentage explained we'd expect from the literature

#interesting that excluding couchii helps mutation-based model more
# than allele-based model

#6: Pulling out specific data for chris ----
#I realised all the stuff above was done with a combined multi-species dataset,
# so i need to pull out just the couchii individuals
couchiidat <- multidat[multidat$species == "Thamnophis couchii",]
summary(couchiidat)
#IDs of individuals we don't have masses for (but do have MAMU for)
nomass <- couchiidat[is.na(couchiidat$mass_g),]
nomass
#I was wrong and we actually do have masses for all of them, which is great

#IDs of individuals we don't have MAMU for
#realized I have MAMUs for some of these so I need to merge dataframes
oldmamu <- read.csv("old_MAMUs.csv")
oldmamu <- oldmamu[,c(2,3,9)]
oldmamu$MAMU <- as.numeric(oldmamu$MAMU)
#some of this has the same collpreno for unknowable reasons so i'm just keeping
# the first instance of each
collpreno <- oldmamu[!duplicated(oldmamu$CollPreNo), ]
indiv <- oldmamu[!duplicated(oldmamu$INDIV), ]

#update MAMU info
couchiidat <- couchiidat %>% 
  filter(!is.na(CollPreNo)) %>% 
  rows_update(distinct(collpreno), by = "CollPreNo", unmatched = "ignore") %>% 
  rows_update(distinct(indiv), by = "INDIV", unmatched = "ignore")
nomamu <- couchiidat[is.na(couchiidat$MAMU),]
head(nomamu)

nomamu <- nomamu[,c(2,3,4)]
write.csv(nomamu, file="couchiiwithoutresistance.csv")

#write.csv(couchiidat,"couchiionlydata.csv")
# saving this since i've filled in so many MAMUs now

#IDs of the juveniles that are weirdly fast 
#first need to re-plot them just to make sure they're couchii (i think they were)
ggplot(couchiidat, aes(mass_g, baseline_ms, color = species))+
  geom_point() + geom_smooth()
conditions <- (couchiidat$mass_g < 25 & couchiidat$baseline_ms > 1)
fastcouchii <- couchiidat[conditions,]
fastcouchii$INDIV
fastcouchii$CollPreNo
write.csv(fastcouchii, file="fastcouchiineonates.csv")
