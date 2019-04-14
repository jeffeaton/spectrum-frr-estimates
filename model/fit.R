library(tidyverse)
library(rstan)
rstan_options(auto_write=TRUE)

#'  Load processed data
load(here::here("data", "data.rda"))

moddat <- moddat %>%
  mutate(agegr = fct_collapse(agegr, "40-49" = c("40-44", "45-49")),
         agegrhiv = fct_collapse(agegr, "35-49" = c("35-39", "40-49")),
         age15to19hiv = as.integer((agegr == "15-19") * hivstatus),
         country = factor(country))

source(here::here("model", "functions.R"))


#' Compile Stan model

sm <- stan_model('model.stan')

#' Prepare data inputs

sdat_all <- make_sdat(dat = moddat)
sdat_bef15 <- make_sdat(dat = subset(moddat, survyear <= 2015))
sdat_fert <- make_sdat(f = ~-1 + surveyid:agegr + yearbefore:hivstatus:agegrhiv,
                       dat = subset(moddat, outcome == "asfr"))
sdat_currpreg <- make_sdat(f = ~-1 + surveyid:agegr,
                           dat = subset(moddat, outcome == "currpreg"))

#' Fit model
fit_all <- sdat_all %>%
  list(dat = .,
       fit = sampling(sm, data=., chains=12, cores=12, iter=1000, refresh=5))

fit_bef15 <- sdat_bef15 %>%
  list(dat = .,
       fit = sampling(sm, data=., chains=12, cores=12, iter=1000, refresh=5))

fit_fert <- sdat_fert %>%
  list(dat = .,
       fit = sampling(sm, data=., chains=12, cores=12, iter=1000, refresh=5))

fit_currpreg <- sdat_currpreg %>%
  list(dat = .,
       fit = sampling(sm, data=., chains=12, cores=12, iter=1000, refresh=5))


#' Save model fits

saveRDS(fit_all, here::here("model", "fits", "fit_all.rds"))
saveRDS(fit_bef15, here::here("model", "fits", "fit_bef15.rds"))
saveRDS(fit_fert, here::here("model", "fits", "fit_fert.rds"))
saveRDS(fit_currpreg, here::here("model", "fits", "fit_currpreg.rds"))
