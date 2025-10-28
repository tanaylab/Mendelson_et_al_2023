# ---
# jupyter:
#   jupytext:
#     formats: ipynb,Rmd,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R 4.0.3
#     language: R
#     name: ir
# ---

# + [markdown] tags=[]
# # Building age-dependent longevity models
# In this notebook we will describe the method for training a probablistic model for reaching the age of 85.  
# The code below exemplifies the mldpEHR package which is used for training age-dependent multivariate prediciton models for longevity and disease.
# The examples will use datasets available in the mldpEHR.data package.
#                                     
#                                                 
#
# -

# ## Setting up simulation data
# The procedures demonstrated in this notebook will used a simulated dataset of patients stored mldpEHR.data dataset.  
# The longevity simulated data consists of the following:
# * longevity.patients - a list of data frames, one for each age, containing the entire population of patients
# * longevity.features - a list of data frames, one for each age, containing the features to be used for training the prediction models.

# +
# installing mldpEHR, and mdlrEHR.data
#remotes::install_github("tanaylab/mldpEHR")
#remotes::install_github("tanaylab/mldpEHR.data")
# -

library(tidyverse)
library(mldpEHR)
library(mldpEHR.data)
names(longevity.patients)
head(longevity.patients[["80"]])

# ### defining model parameters

SURVIVAL_YEARS <- 5
STEP <- 5
MAX_MISSING_PER_FEATURE <- 0.8
FOLDS <- 5

# + [markdown] tags=[]
# ## Building longevity models
# For each age, we will build a classification model for patient survival to age>= 85. 
# As our followup time is limited, we will use older age model score to define the target classification for the younger age model, basically stitching these models together.
# The data we need to provide consists of the entire popultation at each given age along with their sex, age, age at death and potential followup time (time until the end of the database).
#
#
#
# -

longevity <- mldpEHR.mortality_multi_age_predictors(longevity.patients, longevity.features, step=5, nfolds=5, required_conditions='has_cbc')

# + [markdown] tags=[]
# ### Looking at feature significance
#
#
# -

features_sig <- purrr::map(longevity, ~ mldpEHR.prediction_model_features(.x)$summary %>% arrange(desc(mean_abs_shap)))
head(features_sig[[1]])

N_PATIENTS <- 10000
shap_features_80 <- mldpEHR.prediction_model_features(longevity[["80"]])$shap_by_patient %>% 
                             filter(feature %in% head(features_sig[["80"]] %>% pull(feature))) %>% 
                             group_by(feature) %>% 
                             sample_n(N_PATIENTS) %>% 
                             ungroup %>% 
                             mutate(feature=factor(feature, levels=head(features_sig[["80"]] %>% pull(feature))))
options(repr.plot.width=14, repr.plot.height=2.5)
ggplot(shap_features_80, aes(x=value, y=shap)) + geom_point(size=0.01, alpha=0.3) + facet_wrap(~feature, nrow=1, scales="free_y") + theme_bw()

# ## Computing Markovian probability model
#

longevity_markov <- mldpEHR.mortality_markov(longevity, SURVIVAL_YEARS, STEP, seq(0, 1, by=0.1), required_conditions=glue::glue("time >= as.Date('2005-01-01') & time < as.Date('2016-01-01')"))


longevity_prob <- purrr::map2_df(longevity_markov, names(longevity_markov), ~ as_tibble(.x$model[[1]], rownames='sbin') %>% mutate(sex='male', age=.y) %>% bind_rows(as_tibble(.x$model[[2]], rownames='sbin') %>% mutate(sex='female', age=.y)))
options(repr.plot.width=14, repr.plot.height=2.5)
ggplot(longevity_prob %>% mutate(sbin=factor(sbin, levels=c(1:10, "death", "no_score"))), 
       aes(x=sbin, y=death, colour=factor(sex))) + geom_point() + facet_grid(.~age) + theme_bw()


# # Build a disease model for diabetes
# similar to longevity, will used simulated diabetes data, found in mldpEHR.data dataset:
# * diabetes.patients - a list of data frames, one for each age, containing the entire population of patients
# * diabetes.features - a list of data frames, one for each age, containing the features to be used for training the prediction models.
#

diabetes <- mldpEHR.disease_multi_age_predictors(diabetes.patients, diabetes.features, step=5, nfolds=5, required_conditions='has_cbc')

# ##Looking at feature significance

features_sig <- purrr::map2(diabetes, names(diabetes), ~ mldpEHR.prediction_model_features(.x)$summary %>% mutate(age=.y) %>% 
                            arrange(desc(mean_abs_shap)))
head(features_sig[[1]])


N_PATIENTS <- 10000
shap_features_50 <- mldpEHR.prediction_model_features(diabetes[["50"]])$shap_by_patient %>% 
                             filter(feature %in% head(features_sig[["50"]] %>% pull(feature))) %>% 
                             group_by(feature) %>% 
                             sample_n(N_PATIENTS) %>% 
                             ungroup %>% 
                             mutate(feature=factor(feature, levels=head(features_sig[["50"]] %>% pull(feature))))
options(repr.plot.width=14, repr.plot.height=2.5)
ggplot(shap_features_50, aes(x=value, y=shap)) + geom_point(size=0.01, alpha=0.3) + facet_wrap(~feature, nrow=1, scales="free_y") + theme_bw()

# +
## Computing Markovian probability model
# -

diabetes_markov <- mldpEHR.disease_markov(diabetes, 5, 5, seq(0, 1, by=0.1), required_conditions=glue::glue("time >= as.Date('2005-01-01') & time < as.Date('2016-01-01')"))

diabetes_prob <- purrr::map2_df(diabetes_markov, names(diabetes_markov), ~ 
    as_tibble(.x$model[[1]], rownames='sbin') %>% mutate(sex='male', age=.y) %>% 
    bind_rows(
        as_tibble(.x$model[[2]], rownames='sbin') %>% mutate(sex='female', age=.y))
    ) %>% 
    mutate(sbin=factor(sbin, levels=c(1:10, "disease", "disease_death", "death", "no_score")),
          total_disease=disease+disease_death)
options(repr.plot.width=14, repr.plot.height=2.5)
ggplot(diabetes_prob %>% filter(as.numeric(sbin) <= 10), aes(x=sbin, y=total_disease, colour=factor(sex))) + geom_point() + facet_grid(.~age) + theme_bw()

shap <- mldpEHR.prediction_model_features(diabetes[["50"]])

head(shap$shap_by_fold)

nrow(shap$shap_by_fold)


