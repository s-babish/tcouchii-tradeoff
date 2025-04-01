###################################
### Sage Babbish FOCAL analysis ###
###      Due date:: 04-16.      ###
###################################
library(tidyverse)
library(readxl)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load in data
meta <- read.csv("Metadata for multispeciesdata.csv")
data <- read.csv("multispeciesdata.csv")

# a quick summary
summary(data)

##########################
### Data visualization ###
##########################
# Lets do some Quick correlation plots to see how correlated SVL, mass, baseline_ms, sex, species and MAMU are
# library(corrplot)
# corr <- data %>% 
#   drop_na( ) %>% 
#   select(., baseline_ms, MAMU)
# head(corr)
# corrplot <- corrplot(corr = corr, method = "square", type = "upper")
### OMG there are SOOOO many NAs hahaha

# goping to try and Impute
#install.packages("mice")
library(mice)
data_numeric <- data %>% 
  select(., mass_g, SVL_cm, TL_cm, baseline_ms, MAMU)

md.pattern(data_numeric) # this will produce a visual representation of the missing values we are going to need to impute

data_numeric <- data_numeric %>% #### After looking at the TL_cm data, most of it is missing do I'm not going to impute it 
  select(., mass_g, SVL_cm, MAMU, baseline_ms ) %>% 
  rename(., MASS = mass_g, SVL = SVL_cm, SPRINT = baseline_ms)
head(data_numeric)

data_imputed <- data.frame(
  original = data_numeric$SVL,
  imputed_pmm = complete(mice(data_numeric, method = "pmm"))$SVL,
  imputed_cart = complete(mice(data_numeric, method = "cart"))$SVL,
  imputed_lasso = complete(mice(data_numeric, method = "lasso.norm"))$SVL
) #### this will produce a warning of "logged events" which i'm going to ignore.
# The short explantion of the "loged events warning is 
head(data_imputed)

par(mfrow=c(2,2))
hist(data_imputed$original)
hist(data_imputed$imputed_pmm, col = "blue")
hist(data_imputed$imputed_cart, col = "red") # this looks like the closest to the OG data to me
hist(data_imputed$imputed_lasso, col = "green")
par(mfrow = c(1,1))
# Okay looks good, and we will use "cart" (classification and regression trees) method 


# set up all data as numeric
imp.meth <- c("","cart","cart", "" )
data_numeric <- mutate(data_numeric, MAMU = as.numeric(MAMU))

data_imputed <- mice(data = data_numeric, method = imp.meth, 
                     m=30 , maxit = 5, seed = 1234)
summary(data_imputed)

plot(data_imputed) ### Looks good, no notable trends here. 
# it would be bad to see the mean or sd trending up or down consistently 

# next we need to check to see if the imputed data look like reasonable values compared to the real data
stripplot(data_imputed) # blue dots are real data, red dots are imputed values
# also looks good! On to the SEM



#################################
###           SEM             ###
#################################

#First we need a covariance matrix
df.combined <- data_numeric
# Get the imputed data in a workable format
imp.SVL <- apply(data_imputed$imp$SVL, MARGIN = 1, FUN = mean) # IDK if this is right but I'm running out of time and i figure it's better to walk through
#the rest of the code for the SEM and to have a model run, than worry about how to deal with collapsing imputed data
imp.MAMU <- apply(data_imputed$imp$MAMU, MARGIN = 1, FUN = mean)

df.combined$SVL[is.na(df.combined$SVL)] <- imp.SVL # i'm sure there is a better way to do this, but replace NA's with imputed values 
df.combined$MAMU[is.na(df.combined$MAMU)] <- imp.MAMU

#re-check for NA's
which(is.na(df.combined)) # NO NA's LETS GOOOOO

cov.mat <- cov(df.combined)
print(cov.mat) # well I think that this is messed up, and since we imputed the data, this would make sense? 
#But IDK why the diagonals are not 100.... 




#For reference, levaan syntax!!!!

# ~ predict, used for regression of observed outcome to observed predictors (e.g., y ~ x)
# =~ indicator, used for latent variable to observed indicator in factor analysis measurement models (e.g., f =~ q + r + s)
# ~~ covariance (e.g., x ~~ x)
# ~1 intercept or mean (e.g., x ~ 1 estimates the mean of variable x)
# 1* fixes parameter or loading to one (e.g., f =~ 1*q)
# NA* frees parameter or loading (useful to override default marker method, (e.g., f =~ NA*q)
# a* labels the parameter ‘a’, used for model constraints (e.g., f =~ a*q)
                                
#In this model, each of the variables (MAMU, SVL, and MASS) directly influences the response variable SPRINT.

#install.packages(c("levaan","semPlot"))
library(lavaan)

# Define the model
model1 <- '
  
  # Structural model
  SPRINT ~ MAMU
  SPRINT ~ SVL
  SPRINT ~ MASS
  
'

# Fit the model to the data
SEM1 <- sem(model1, data = df.combined)

# Summarize the results
summary(SEM1, fit.measures = TRUE)

# Load the semTools package
#install.packages("lavaanPlot")
library(lavaanPlot)

# Create a path diagram for model 1
lavaanPlot(SEM1)


###### Well thats just not great...

# IDK what to do now, so hopefully this was somewhat helpful. Apologies if it is useless. 
