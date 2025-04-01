##### Analysis of Sage Babish's Data ##### 
##### Author: Victoria Peechatt 
##### Date: 04/17/2024

##### Objectives for Sage: #####
##### -  Figure out direct and indirect associations between variables 
##### -  relationships between alleles, SVL, MAMU, and sprint speed 
##### -  Modeling predictors of TTX resistance

##### Objectives for myself: #####
##### - Learn how to SEM
##### - and test different causal path models  
##### - Learn about data imputation 

##### Package Test ######

#checks to see if you have package
#if not, installs it 

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

pkgTest("data.table")
pkgTest("inspectdf")
pkgTest("tidyverse")
pkgTest("gridExtra")
pkgTest("lme4")
pkgTest("sjPlot")
pkgTest("effects")
pkgTest("vegan")
pkgTest("asbio")
pkgTest("lavaan")
pkgTest("lavaanPlot")
pkgTest("foreign")
pkgTest("mice")

##### Loading + Visualizing Data #####

babish.data = data.table::fread(
    "sbabish_data.csv",
)
head(babish.data)

type_info = inspectdf::inspect_types(babish.data)
type_info$col_name
df = inspectdf::inspect_cat(babish.data)
inspectdf::show_plot(df)
inspectdf::show_plot(df, high_cardinality=1)

#####
# Making sure MAMU is numeric (there's notes in some of the cells)
# Also asked Sage about the large values of MAMU
# and they said that the large values have to be data entry errors 
# so filtering out MAMU over 500, change line below to keep desired values 
#####

babish.data$MAMU = as.numeric(babish.data$MAMU)
babish.filtered = babish.data %>% filter(babish.data$MAMU < 500)
babish.filtered$age = as.factor(babish.filtered$age)
babish.filtered$sex = as.factor(babish.filtered$sex)
babish.filtered$allele_a = as.factor(babish.filtered$allele_a)
babish.filtered$allele_b = as.factor(babish.filtered$allele_b)
babish.filtered$species = as.factor(babish.filtered$species)

str(babish.filtered)

##### Imputation #####
# Imputing data using logreg for sex which uses logistic regression 
# and norm for SVL_cm and TL_cm which uses bayesian linear regression 
# and pmm for all the others which uses predictive mean matching
# see ?mice or methods(mice) for more details 
#####

# Shows how much of each column is missing 
data_missing <- unlist(lapply(babish.filtered, function(x) sum(is.na(x))))/nrow(dataset)
sort(data_missing[data_missing > 0], decreasing = TRUE)

# mice mice baby 
imputation1 <- mice(babish.filtered, method="pmm", maxit=5)
meth = imputation1$method
meth["SVL_cm"] = "norm"
meth["TL_cm"] = "norm"
meth["sex"] = "logreg"
imputation1 <- mice(dataset, method=meth, maxit=5)

# Shows what variables are being used to impute others 
# You can change the predictors by changing the pred object 

imputation1$predictorMatrix

# Inspect the convergence of the mice algorithm
# Looks like they're all overlapping and no consistent trend in any chain 
# so thats good 

plot(imputation1)

# The imputed data: 

#imputation1$imp
head(imputation1$imp$SVL_cm)

# Reincorporating it into the rest of the data, if you don't want to use 
# the imputed data delete the following line

babish.filtered = complete(imputation1, 1)

##### Plotting as an initial look

ggplot(babish.filtered, aes(MAMU, baseline_ms))+
  geom_point()+
  geom_smooth(method = lm)

#####
# Looks like a log normal transformation of MAMU would do good things, thanks Lee 
# Transformation of MAMU 1 - log, 2 - standardize 
######

babish.filtered$MAMU_log = log(babish.filtered$MAMU)
babish.filtered$mass_g = log(babish.filtered$mass_g)
babish.filtered$SVL_cm = log(babish.filtered$SVL_cm)

ggplot(babish.filtered, aes(MAMU_log, baseline_ms))+
  geom_point() + 
  geom_smooth()

# Standardizing MAMU manually 
# doing it in one line would look like this: 
# variables_to_std <- decostand(variables_to_std, method = "standardize")
babish.filtered$MAMU_stand = (babish.filtered$MAMU_log - 
                                mean(babish.filtered$MAMU_log, na.rm = TRUE))/
  sd(babish.filtered$MAMU_log, na.rm = TRUE)

# Looking at more plots 

ggplot(babish.filtered, 
       aes(SVL_cm, baseline_ms, 
           color = species))+
  geom_point()+
  geom_smooth(method = lm)

ggplot(babish.filtered, 
       aes(MAMU_stand, baseline_ms, 
           color = allele_a))+
  geom_point()+
  geom_smooth(method = lm)

#####
# Species and alleles seem to group together in plots 
# so making them dummy variables next 
##### 
babish.filtered.split = babish.filtered %>% mutate(species = str_split(species, " "))
babish.filtered.split$species = babish.filtered.split$species %>% lapply(function(x) x[2])

babish.filtered.split = mutate(babish.filtered.split, ones = 1) |> 
  spread(key = species, value = ones, fill = 0)

babish.filtered.split = mutate(babish.filtered.split, ones = 1) |> 
  spread(key = allele_a, value = ones, fill = 0, sep = "_")

babish.filtered.split = mutate(babish.filtered.split, ones = 1) |> 
  spread(key = allele_b, value = ones, fill = 0, sep = "_")

# and standardizing the variables 
# keeping split dummy variables data as separate 
data_scaled = apply(babish.filtered[,c(9:12,16)], MARGIN = 2, scale)
data_scaled = cbind(data_scaled, babish.filtered[,c(1,7,8,14,15)])

data_split_scaled = apply(babish.filtered.split[,c(8:13)], MARGIN = 2, scale)
data_split_scaled = cbind(data_split_scaled, babish.filtered.split[,c(15:34)])

##### Generalized Linear mixed regression models ##### 

# Model 1 

model1 = lmer(baseline_ms ~ MAMU_log + mass_g*SVL_cm + (1|species), data = data_scaled)

summary(model1)

# Model 2 

model2 = lmer(baseline_ms ~ MAMU_log + age + mass_g + SVL_cm + TL_cm + 
                SVL_cm*TL_cm + SVL_cm*mass_g + 
                age*SVL_cm + species + (1|sex), data = data_scaled)

summary(model2)

plot_model(model2, 
           show.values = TRUE, show.p = FALSE,
           title = "Effect of Variables on Baseline Speed")
tab_model(model2, show.re.var = TRUE, show.p = FALSE)

effect_plot = effect(term = "MAMU_log", model2)
effect_plot_df = as.data.frame(effect_plot)

model2_plot = ggplot() + 
  geom_point(data = data_scaled, aes(MAMU_log, baseline_ms, color = species)) + 
  geom_point(data = effect_plot_df, aes(x = MAMU_log, y= fit), color = "blue") + 
  geom_line(data = effect_plot_df, aes(x = MAMU_log, y= fit), color = "blue") + 
  geom_ribbon(data = effect_plot_df, 
              aes(x = MAMU_log, ymin = lower, ymax = upper),
              alpha = 0.3,
              fill = "blue")
model2_plot


##### SEM's woooooooooo ##### 

# Model 3

model3 = '
  mass_g ~ SVL_cm
  baseline_ms ~ MAMU_log + SVL_cm + mass_g'

model3.fit = sem(model3, data = data_scaled)
lavResiduals(model3.fit)
summary(model3.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "model3", model = model3.fit, coefs = TRUE)


##### Model 3 Interpretation #####
# RMSEA: 0.000
# Chi-square: 0.594
# df: 1
# This is a shnazzy fit 
# Interpretation: 
# SVL has a pretty much direct, positive standardized effect on mass 
# which then has a positive 0.31 standardized effect on speed 
# MAMU_log has a 0.16 standardized effect on speed

# Model 4

model4 = '
  baseline_ms ~ MAMU_log + SVL_cm + mass_g + species
  MAMU_log ~ SVL_cm + mass_g + species
  mass_g ~ sex + species'

# fitting model to data 

model4.fit = sem(model4, data = data_scaled, cluster = "species")

summary(model4.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "model4", model = model4.fit, coefs = TRUE)


##### Model 4 Interpretation #####
# RMSEA: 0.425
# Chi-square: 415.683
# df: 3
# This is an aight fit 
# Interpretation: 
# Key takeaways are that sex has an effect on mass which has an effect on speed
# everything else looks meh 
# but also it's hard to interpret with sex and species being categorical 
# *foreshadowing the future of this code...*


# Model 5

model5 = '
  mass_g ~ SVL_cm
  MAMU_log ~ species
  baseline_ms ~ MAMU_log + SVL_cm + mass_g'

model5.fit = sem(model5, data = data_scaled)

summary(model5.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "model5", model = model5.fit, coefs = TRUE)

##### Model 5 Interpretation #####
# RMSEA: 0.127
# Chi-square: 52.941
# df: 4
# This is an excellent fit 
# Interpretation: 
# SVL has a positive effect on mass (makes sense)
# and a slightly negative effect on baseline speed 
# MAMU_log has a positive standardized effect on speed (0.16) 
# while mass has double that effect on speed
# species seems to be effecting MAMU_log but unsure what it means since its a fac
# good thing we split it up for lter 


##### 
# Using split up dummy variable data 
# Can't add all the categories to a single model because it should be k-1 groups 
#####

model6 = '
  MAMU_log ~ elegans + atratus + sirtalis
  baseline_ms ~ MAMU_log + elegans + atratus + sirtalis
'

model6.fit = sem(model6, data = data_split_scaled)

summary(model6.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "model6", model = model6.fit, coefs = TRUE)

##### Model 6 Interpretation #####
# RMSEA: 0
# Chi-square: 0
# df: 0
# This is a perfect fit, oversaturated  
# Interpretation: The elegans and atratus species have a negative standardized 
# effect (compared to couchii) on MAMU_log while sirtalis has a slightly 
# more positive effect on MAMU_log. All three species have a negative std 
# effect on speed compared to couchii 

# Model 7 

model7 = '
  MAMU_log ~ allele_a_EPN + allele_a_LVNV + allele_a_P + allele_a_T + 
  allele_a_V +  allele_a_WT
  baseline_ms ~ MAMU_log
'
model7.fit = sem(model7, data = data_split_scaled)

summary(model7.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "model7", model = model7.fit, coefs = TRUE)


##### Model 7 Interpretation #####
# RMSEA: 0.221
# Chi-square: 228.624
# df: 6
# This is a pretty good fit 
# Interpretation: The T, V, P, and WT alleles have a negative standardized 
# effect (compared to an "NA" for allele A) on MAMU_log
# THE EPN and LVNV alleles have a positive standardized effect on MAMU_log 


model8 = '
  MAMU_log ~ allele_a_EPN + allele_a_LVNV + allele_a_P + allele_a_T + 
  allele_a_V +  allele_a_WT
  baseline_ms ~ allele_a_EPN + allele_a_LVNV + allele_a_P + allele_a_T + 
  allele_a_V +  allele_a_WT
  baseline_ms ~ MAMU_log
'
model8.fit = sem(model8, data = data_split_scaled)

summary(model8.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE)

lavaanPlot(name = "model8", model = model8.fit, coefs = TRUE)

##### Model 8 Interpretation #####
# RMSEA: 0.0000
# Chi-square: 0.0000
# df: 0
# This is a perfectly fit, saturated model 
# Interpretation: The T allele is directly relational to baseline speed, 
# having a positive standardized effect (compared to an "NA" for allele A), 
# while having a negative standardized effect of -0.63 on MAMU_log
# Having P or WT as allele A has a strong negative standardized  effect on MAMU 
# (-1.45 and -1.34 respectively) on MAMU_log while having a positive std effect 
# on baseline speed (0.27 and 0.24 respectively)
# The other alleles seem to have an influence on MAMU_log but not baseline speed
