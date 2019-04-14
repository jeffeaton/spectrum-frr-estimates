library(tidyverse)
library(rdhs)
library(haven)
library(here)
library(survey)
library(demogsurv) # devtools::install_github("mrc-ide/demogsurv@2e9b4f0")

#' Identify country codes for desired countries
cc <- rdhs::dhs_countries() %>%
  filter(RegionName == "Sub-Saharan Africa") %>%
  .$DHS_CountryCode

#' SurveyCharacteristicID = 23 is HIV testing
dhs_survey_characteristics()

#' Identify all surveys from SSA
surveys <- dhs_surveys(countryIds = cc, surveyYearStart=2000) %>%
  filter(RegionName == "Sub-Saharan Africa")
  

#' Identify datasets for these surveys
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
ard <- dhs_datasets(fileType = "AR", fileFormat = "flat", surveyIds = surveys$SurveyId)

#' Get local path to dataset (download if needed)
ird$path <- unlist(get_datasets(ird))
ard$path <- unlist(get_datasets(ard))

#' Variables to extract from individual recode datasets
#' * v213: currently pregnant
#' * v525: age at first sex
#' * v527: time since last sex
#' 

ir <- readRDS(tail(ird$path, 1))

bvars <- expand.grid("b",
                     c("idx", "ord", 1:7),
                     "_",
                     formatC(1:20, digits=1, format="d", flag="0")) %>%
  {do.call(mapply, c(paste0, .))}

#' 
irvars <- c("caseid", "v001", "v002", "v003", "v005", "v008", "aidsex", "v011", "v012",
            "v021", "v022", "v023", "v024", "v025", "v213", "v525", "v527", bvars)
arvars <- c("hivclust", "hivnumb", "hivline", "hiv03", "hiv05", "shiv51")
allvars <- c("surveyid", "country", "survyear", irvars, "hiv03", "hiv05", "shiv51")

            
#' Load and merge datasets
datlst <- list()
for(survid in ird$SurveyId){

  print(survid)
  
  ir <- filter(ird, SurveyId == survid)$path %>% readRDS

  if(survid == "CI2005AIS"){
    ## For Cote d'Ivoire 2005 AIS, individuals are uniquely identified by
    ## four variables {cluster, structure, household, line}.
    ir$v002 <- 100*ir$sstruct + ir$v002
  }

  if(!exists("aidsex", ir)){
    ir$aidsex <- haven::labelled(2, c("men" = 1, "women" = 2), label = "Sex")
  }
  
  dat <- ir %>%
    select(intersect(irvars, names(ir))) %>%
    filter(aidsex == 2)
  dat[setdiff(irvars, names(ir))] <- NA
  
  if(nrow(ard[SurveyId == survid])){
    ar <- filter(ard, SurveyId == survid)$path %>% readRDS
    if(survid == "CI2005AIS")
      ar$hivnumb <- 100*ar$hivstruct + ar$hivnumb
    ar[setdiff(arvars, names(ar))] <- NA

    dat <- dat %>%
      left_join(ar %>% select(arvars),
                by=c("v001" = "hivclust",
                     "v002" = "hivnumb",
                     "v003" = "hivline"))
  }

  dat[setdiff(allvars, names(dat))] <- NA
  
  dat$surveyid <- survid
  dat$country <- filter(ird, SurveyId == survid)$CountryName
  dat$survyear <- filter(ird, SurveyId == survid)$SurveyYear

  datlst[[survid]] <- dat %>% select(allvars)
}

## Replace results for Zambia 2013-14 with confirmed test results

datlst[["ZM2013DHS"]]$hiv03 <- datlst[["ZM2013DHS"]]$shiv51

## Create HIV, pregnancy status, and sexual history indicators

datlst <- lapply(datlst, mutate,
                 prev = if_else(hiv03 %in% 0:3, as.integer(hiv03 > 0), NA_integer_),
                 hivstatus = factor(prev, 0:1, c("negative", "positive")),
                 currpreg = if_else(v213 %in% 0:1, v213, NA_integer_),
                 eversex = as.integer(v525 != 0),
                 recentsex = as.integer(eversex & if_else(!v527 %in% 997:999, v527 <= 204 | v527 == 995, NA)),
                 sex12m = as.integer(eversex & if_else(!v527 %in% 997:999, v527 <= 400 | v527 == 995, NA)))


#' Calculate indicators

library(parallel)
options(mc.cores = parallel::detectCores())

## Surveys with birth histories
hasbh <- !names(datlst) %in% c("TZ2003AIS", "UG2011AIS", "CG2009AIS", "MZ2009AIS")


tfr <- mclapply(datlst[hasbh], calc_tfr,
                by=~surveyid + country + survyear,
                tips=0:6, strata=NULL)

asfr <- mclapply(datlst[hasbh], calc_asfr,
                 by=~surveyid + country + survyear,
                 tips=0:6, strata=NULL, counts=TRUE)

gfr <- mclapply(datlst[hasbh], calc_asfr,
                 by=~surveyid + country + survyear,
                 tips=0:6, agegr=c(15, 60), strata=NULL, counts=TRUE)


## ASFR by HIV status
hashiv <- sapply(datlst, function(x) sum(!is.na(x$hivstatus)) > 0)

datlst_hiv <- datlst[hashiv & hasbh] %>%
  lapply(filter, !is.na(hiv05))

tfrhiv <- mclapply(datlst_hiv, calc_tfr,
                   by=~surveyid + country + survyear + hivstatus,
                   weight = "hiv05",
                   tips=0:6, strata=NULL)

asfrhiv <- mclapply(datlst_hiv, calc_asfr,
                    by=~surveyid + country + survyear+hivstatus,
                    weight = "hiv05",
                    tips=0:6, strata=NULL, counts=TRUE)

gfrhiv <- mclapply(datlst_hiv, calc_asfr,
                   by=~surveyid + country + survyear+hivstatus,
                   weight = "hiv05", agegr=c(15, 50),
                   tips=0:6, strata=NULL, counts=TRUE)


## Percentage currently pregnant

calc_bin <- function(dat,
                     agegr=c(15, 50),
                     byhiv=FALSE,
                     weights=if(byhiv) ~hiv05 else ~v005,
                     strata=NULL,
                     formula=~currpreg){

  dat$agegr <- cut(dat$v012, agegr, demogsurv:::.epis_labels(agegr), TRUE, FALSE)
  if(!byhiv)
    dat$hivstatus <- "all"
  des <- svydesign(ids=~v021, data=dat[!is.na(dat[[all.vars(weights)]]), ],
                   strata=strata, weights=weights)

  val <- left_join(svyby(formula, ~surveyid + country + survyear + hivstatus + agegr,
                         des, unwtd.count, na.rm=TRUE) %>%
                   rename(n = counts) %>% select(-se),
                   svyby(formula, ~surveyid + country + survyear + hivstatus + agegr,
                         des, svyciprop, vartype=c("se", "ci"), na.rm=TRUE))
  names(val) <- sub("^se\\..*", "se", names(val))
  
  val
}

dhs_preg <- bind_rows(
  mclapply(datlst, calc_bin),
  mclapply(datlst, calc_bin, agegr=3:10*5),
  mclapply(datlst[hashiv], calc_bin, byhiv=TRUE),
  mclapply(datlst[hashiv], calc_bin, agegr=3:10*5, byhiv=TRUE)
)

## Ever had sex

hasevsx <- sapply(datlst, function(x) sum(!is.na(x$eversex)) > 0)

dhs_eversex <- bind_rows(
  mclapply(datlst[hasevsx], calc_bin, formula=~eversex),
  mclapply(datlst[hasevsx], calc_bin, agegr=3:10*5, formula=~eversex),
  mclapply(datlst[hasevsx & hashiv], calc_bin, byhiv=TRUE, formula=~eversex),
  mclapply(datlst[hasevsx & hashiv], calc_bin, agegr=3:10*5, byhiv=TRUE, formula=~eversex)
)

dhs_sex12m <- bind_rows(
  mclapply(datlst[hasevsx], calc_bin, formula=~sex12m),
  mclapply(datlst[hasevsx], calc_bin, agegr=3:10*5, formula=~sex12m),
  mclapply(datlst[hasevsx & hashiv], calc_bin, byhiv=TRUE, formula=~sex12m),
  mclapply(datlst[hasevsx & hashiv], calc_bin, agegr=3:10*5, byhiv=TRUE, formula=~sex12m)
)


## Prevalence among currently pregnant

datpreg <- lapply(datlst[hashiv], subset, currpreg & v012 %in% 15:49)

dhs_pregprev <- bind_rows(
  mclapply(datpreg, calc_bin, formula=~prev, weights=~hiv05),
  mclapply(datpreg, calc_bin, agegr=3:10*5, formula=~prev, weights=~hiv05)
) %>%
  mutate(hivstatus = NULL)

dhs_prev <- bind_rows(
  mclapply(datlst[hashiv], calc_bin, formula=~prev, weights=~hiv05),
  mclapply(datlst[hashiv], calc_bin, agegr=3:10*5, formula=~prev, weights=~hiv05)
) %>%
  mutate(hivstatus = NULL)


## pool and save

dhs_tfr <- bind_rows(tfrhiv, tfr %>% bind_rows %>% mutate(hivstatus = "all")) %>%
  rename(se = se_tfr)
dhs_asfr <- bind_rows(asfrhiv, asfr %>% bind_rows %>% mutate(hivstatus = "all")) %>%
  rename(se = se_asfr)
dhs_gfr <- bind_rows(gfrhiv, gfr %>% bind_rows %>% mutate(hivstatus = "all")) %>%
  rename(se = se_asfr)

save(dhs_preg, dhs_eversex, dhs_sex12m, dhs_pregprev, dhs_prev, dhs_tfr, dhs_asfr, dhs_gfr,
     file = here("data", "dhs.rda"))
