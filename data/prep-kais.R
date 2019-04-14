library(tidyverse)
library(haven)
library(lubridate)
library(survey)
library(demogsurv)



#' Recode KAIS 2012 data

ir <- as.data.frame(read_dta("~/Documents/Data/Kenya/KAIS/KAIS 2012/alladults.dta"))
br <- as.data.frame(read_dta("~/Documents/Data/Kenya/KAIS/KAIS 2012/allbirths.dta"))

ir <- ir %>%
  mutate(intvd = decimal_date(as.Date(interviewdate, origin = "1960-01-01")),
         intvc = floor(12 * (intvd - 1900)),
         age = q_102,
         agegr = cut(age, 3:10*5, paste0(3:9*5, "-", 3:9*5+4), TRUE, FALSE),
         sex = factor(sex, 1:2, c("male", "female")),
         dob_d = as.integer(sub("(.+)/(.+)/(.+)", "\\1", q_101)),
         dob_m = as.integer(sub("(.+)/(.+)/(.+)", "\\2", q_101)),
         dob_y = as.integer(sub("(.+)/(.+)/(.+)", "\\3", q_101)))

## Impute missing DOB based on age
ir$dob_y[ir$dob_y %in% c(8888, 9898)] <- floor(ir$intvd[ir$dob_y %in% c(8888, 9898)] - (ir$q_102[ir$dob_y %in% c(8888, 9898)] + 0.5))
ir$dob_d <- as.integer(sub("[89]8", "15", ir$dob_d))
ir$dob_m <- as.integer(sub("[89]8", "6", ir$dob_m))

ir <- ir %>% mutate(dobc = 12 * (dob_y - 1900) + dob_m,  # respondent dob (CMC)
                    hivstatus = tolower(nhrl_hivresults),
                    prev = as.integer(hivstatus == "positive"),
                    currpreg = as.integer(q_345 == 1),
                    eversex = if_else(eversex %in% 1:2, as.integer(eversex == 1), NA_integer_),
                    sex12m = if_else(sexactive12mos %in% 1:2, as.integer(sexactive12mos == 1), NA_integer_))

br$chsex <- factor(br$q_302, 1:2, c("male", "female"))
br$chdob <- floor(12 * (decimal_date(as.Date(br$q_303, "%d/%m/%Y")) - 1900))


## Drop records with missing date of interview
## Using HIV weights (abweight) for all analysis
ir <- filter(ir, sex == "female" & !is.na(interviewdate) & !is.na(abweight))

## Drop records from birth history that are not in individual dataset
sum(duplicated(ir$femaleid))
table(!br$femaleid %in% ir$femaleid)
br <- filter(br, femaleid %in% ir$femaleid)


## ASFR and TFR calculation
calc_asfr(ir, bhdata=br, clusters = ~qclust, strata = ~strata2, tips=c(0, 1),
          id = "femaleid", dob = "dobc", intv = "intvc", weight = "abweight",
          bvars = "chdob")

calc_tfr(ir, bhdata=br, clusters = ~qclust, strata = ~strata2, tips=c(0, 1),
         id = "femaleid", dob = "dobc", intv = "intvc", weight = "abweight",
         bvars = "chdob")

kais12_tfr <- bind_rows(
  calc_tfr(mutate(ir, hivstatus="all"), bhdata=br,
           by=~hivstatus, tips=c(0, 1),
           clusters = ~qclust, strata = ~strata2,
           id = "femaleid", dob = "dobc", intv = "intvc",
           weight = "abweight", bvars = "chdob"),
  calc_tfr(ir, bhdata=br, by=~hivstatus, tips=c(0, 1),
           clusters = ~qclust, strata = ~strata2,
           id = "femaleid", dob = "dobc", intv = "intvc",
           weight = "abweight", bvars = "chdob")
)

kais12_asfr <- bind_rows(
  calc_asfr(mutate(ir, hivstatus="all"), bhdata=br,
            by=~hivstatus, tips=c(0, 1),
            clusters = ~qclust, strata = ~strata2,
            id = "femaleid", dob = "dobc", intv = "intvc",
            weight = "abweight", bvars = "chdob", counts=TRUE),
  calc_asfr(ir, bhdata=br,
            by=~hivstatus, tips=c(0, 1),
            clusters = ~qclust, strata = ~strata2,
            id = "femaleid", dob = "dobc", intv = "intvc",
            weight = "abweight", bvars = "chdob", counts=TRUE)
)


## Currently pregnant

des <- svydesign(ids=~qclust, data=ir, strata=~strata2, weights=~abweight)

svyprop <- function(formula, by, design){
  val <- left_join(svyby(formula, by, design, unwtd.count, na.rm=TRUE) %>%
                   rename(n = counts) %>% select(-se),
                   svyby(formula, by, design, svyciprop, vartype=c("se", "ci"), na.rm=TRUE))
  names(val) <- sub("^se\\..*", "se", names(val))
  val
}
  
kais12_preg <- bind_rows(
  svyprop(~currpreg, ~hivstatus + agegr,
          update(subset(des, age %in% 15:49), hivstatus = "all", agegr="15-49")),
  svyprop(~currpreg, ~hivstatus+agegr, update(des, hivstatus="all")),
  svyprop(~currpreg, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), agegr = "15-49")),
  svyprop(~currpreg, ~hivstatus+agegr, des)
)
                           
kais12_eversex <- bind_rows(
  svyprop(~eversex, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), hivstatus = "all", agegr="15-49")),
  svyprop(~eversex, ~hivstatus+agegr, update(des, hivstatus="all")),
  svyprop(~eversex, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), agegr = "15-49")),
  svyprop(~eversex, ~hivstatus+agegr, des)
)

kais12_sex12m <- bind_rows(
  svyprop(~sex12m, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49),
                 hivstatus = "all", agegr="15-49")),
  svyprop(~sex12m, ~hivstatus+agegr,
          update(des, hivstatus="all")),
  svyprop(~sex12m, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), agegr = "15-49")),
  svyprop(~sex12m, ~hivstatus+agegr, des)
)

kais12_pregprev <- bind_rows(
  svyprop(~prev, ~agegr, subset(des, currpreg == 1)),
  svyprop(~prev, ~agegr, update(subset(des, currpreg == 1), agegr = "15-49"))
)

kais12_prev <- bind_rows(
  svyprop(~prev, ~agegr, des),
  svyprop(~prev, ~agegr, update(des, agegr = "15-49"))
)


## Save results

idinfo <- data.frame(surveyid = "KE2012KAIS",
                     country = "Kenya",
                     survyear = 2012)

kais12_preg <- data.frame(idinfo, kais12_preg)
kais12_eversex <- data.frame(idinfo, kais12_eversex)
kais12_sex12m <- data.frame(idinfo, kais12_sex12m)
kais12_asfr <- data.frame(idinfo, kais12_asfr)
kais12_tfr <- data.frame(idinfo, kais12_tfr)
kais12_pregprev <- data.frame(idinfo, kais12_pregprev)
kais12_prev <- data.frame(idinfo, kais12_prev)



#' ## KAIS 2007
#'
#' Includes currently pregnant only, not fertility history
#' q214: "Are you pregnant now"
#' q321: "Age at first intercourse" 
#' 
#' q329u: "When was the last time you had sexual intercourse(Units)" 
#' q329n: "When was the last time you had sexual intercourse (Number)" 
#' q348: "Total number of persons one had had sex with last 12months" 
                                    

kais07 <- foreign::read.spss("~/Documents/Data/Kenya/KAIS/KAIS 2007/Individual Questionnire.sav", to.data.frame = TRUE)
table(kais07$HIV)

with(kais07, table(type.convert(q329n), q329u, useNA="always"))

kais07 <- filter(kais07,
                 qresult == "Completed",
                 sex == "Female",
                 q103 %in% 15:49) %>%
  mutate(country = "Kenya",
         surveyid = "KE2007KAIS",
         survyear = 2007,
         psu = qclust,
         strata = qprov,
         region = qprov,
         indweight = ind_weight,
         hivweight = bl_weight,
         sex = tolower(sex),
         age = type.convert(q103),
         agegr = cut(age, 3:10*5, paste0(3:9*5, "-", 3:9*5+4), TRUE, FALSE),
         hivstatus = c("positive", "negative")[match(HIV, 1:2)],
         prev = as.integer(hivstatus == "positive"),
         evertest = (!is.na(q605) & q605 == "Yes" & q609 == "Yes" |
                     (q605 == "No" | is.na(q605)) & q612 == "Yes" & q616 == "Yes"),
         test12m = evertest & (!is.na(q607) & q607 == "Less than 12months ago" |
                               !is.na(q613) & q613 == "Less than 12months ago"),
         aware = evertest & q705 == "Yes",
         art_self = aware & q707 == "Yes",
         currpreg = c(TRUE, FALSE, FALSE)[match(q214, c("Yes", "No", "Unsure"))],
         eversex = q321 != "Never",
         q329n = type.convert(q329n),
         sex12m = eversex & (q329u == "Days" & q329n <= 365 |
                             q329u == "Weeks" & q329n <= 52 |
                             q329u == "Months" & q329n <= 12 |
                             q329u == "Years" & q329n <= 0))


des <- svydesign(ids=~psu, data=kais07, strata=~strata, weights=~bl_weight)

kais07_preg <- bind_rows(
  svyprop(~currpreg, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), hivstatus = "all", agegr="15-49")),
  svyprop(~currpreg, ~hivstatus+agegr, update(des, hivstatus="all")),
  svyprop(~currpreg, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), agegr = "15-49")),
  svyprop(~currpreg, ~hivstatus+agegr, des)
)

kais07_eversex <- bind_rows(
  svyprop(~eversex, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), hivstatus = "all", agegr="15-49")),
  svyprop(~eversex, ~hivstatus+agegr, update(des, hivstatus="all")),
  svyprop(~eversex, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), agegr = "15-49")),
  svyprop(~eversex, ~hivstatus+agegr, des)
)

kais07_sex12m <- bind_rows(
  svyprop(~sex12m, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), hivstatus = "all", agegr="15-49")),
  svyprop(~sex12m, ~hivstatus+agegr, update(des, hivstatus="all")),
  svyprop(~sex12m, ~hivstatus+agegr,
          update(subset(des, age %in% 15:49), agegr = "15-49")),
  svyprop(~sex12m, ~hivstatus+agegr, des)
)

kais07_pregprev <- bind_rows(
  svyprop(~prev, ~agegr, subset(des, currpreg == 1)),
  svyprop(~prev, ~agegr, update(subset(des, currpreg == 1), agegr = "15-49"))
)

kais07_prev <- bind_rows(
  svyprop(~prev, ~agegr, des),
  svyprop(~prev, ~agegr, update(des, agegr = "15-49"))
)

idinfo <- list(surveyid = "KE2007KAIS",
               country = "Kenya",
               survyear = 2007)

kais07_preg <- data.frame(idinfo, kais07_preg)
kais07_eversex <- data.frame(idinfo, kais07_eversex)
kais07_sex12m <- data.frame(idinfo, kais07_sex12m)
kais07_pregprev <- data.frame(idinfo, kais07_pregprev)
kais07_prev <- data.frame(idinfo, kais07_prev)

kais_preg <- bind_rows(kais07_preg, kais12_preg)
kais_eversex <- bind_rows(kais07_eversex, kais12_eversex)
kais_sex12m <- bind_rows(kais07_sex12m, kais12_sex12m)
kais_asfr <- kais12_asfr %>% rename(se = se_asfr)
kais_tfr <- kais12_tfr %>% rename(se = se_tfr)
kais_pregprev <- bind_rows(kais07_pregprev, kais12_pregprev)
kais_prev <- bind_rows(kais07_prev, kais12_prev)

save(kais_preg, kais_eversex, kais_sex12m, kais_asfr, kais_tfr, kais_pregprev, kais_prev,
     file = here::here("data", "kais.rda"))
