library(tidyverse)
library(haven)
library(survey)


hsrc05 <- read_dta("~/Documents/Data/South Africa/HSRC/2005/adult and youth/SABSSM2005_adult_youth.dta") %>% as.data.frame()

## csa12    "Sexual activity in the last 12 months"

za05 <- hsrc05 %>%
  filter(q1_2 == 2, q1_1 %in% 15:49, finresq == 1) %>%
  transmute(surveyid = "ZA2005HSRC",
            survyear = 2005,
            country = "South Africa",
            ea = ea,
            stratum = province,
            age = q1_1,
            eversex = as.integer(na_if(q5_1, 3) == 1),
            sex12m = as.integer(csa12 == 3),
            currpreg = as.integer(q7_19 == 1),
            hivstatus = factor(hivstat, 1:2, c("positive", "negative")),
            prev = as.integer(hivstatus == "positive"),
            weight = ibreal12,
            hivweight = ibreal1)


## hivstat  "HIV status from the lab"
## q7_18    "Q7_18 Have you ever been pregnant"
## q7_19    "Q7_19 Are you currently pregnant"
## q7_20    "Q7_20 Have you been pregnant in the last 24 months"
## finresfh "Response status to giving BLOOD for HIV test"
## finresq  "Response status to the QUESTIONNAIRE"
## agefive  "Age in five year age brackets"
## marstat  "Marital status"
## ibreal12 "The BENCHMARKED weight for those interviewed irrespective of whether they gave s"
## ibreal1  "The BENCHMARKED weight for those interviewed and also gave specimen"
## eatype   "The EA subtype"
## geotype  "The Geotype of the EA number"
## agec     "The age group of respondent"
## q1_1     "Q1_1 How old were you at your last birthday"
## q5_3     "Q5_3 How old were you when you had sex for the first time"
## q5_4     "Q5_4 How old was the person you had sex with for the first time"
## q1_2     "Q1_2 Sex of the respondent"
## rq28     "RQ28:Have you had sex during the past 12 months?"


hsrc08 <- read_dta("~/Documents/Data/South Africa/HSRC/2008/combined/SABSSM2008_Combined_V02.dta")

za08 <- hsrc08 %>%
  filter(ageyrs %in% 15:49, sex == 2, fresp %in% 1:2) %>%
  transmute(surveyid = "ZA2008HSRC",
            survyear = 2008,
            country = "South Africa",
            ea = ea,
            stratum = prov,
            age = ageyrs,
            eversex = if_else(!rq20 %in% c(0, 3), rq20  == 1, NA) %>% as.integer,
            sex12m = ((eversex == 1 | is.na(eversex))  & (na_if(rq28, 3) == 1)) %>% as.integer,
            everpreg = if_else(!rq61 %in% c(0, 3), rq61  == 1, NA),
            currpreg = if_else(eversex == 0 & (is.na(rq64)), FALSE,
            (everpreg == 1 | is.na(everpreg))  & (rq64 == 1)) %>% as.integer,
            hivstatus = factor(hivstat, 1:0, c("positive", "negative")),
            prev = as.integer(hivstatus == "positive"),
            weight = bcreal12,
            hivweight = bcreal1,
            everpreg = NULL)

hsrc12 <- read_dta("~/Documents/Data/South Africa/HSRC public/2012/SABSSM2012_Adult.dta")

cbind(sapply(hsrc12, attr, "label"))

za12 <- hsrc12 %>%
  filter(Q1_1 %in% 15:49 & Q1_2 == 2 & fresp %in% 1:2) %>%
  transmute(surveyid = "ZA2012HSRC",
            survyear = 2012,
            country = "South Africa",
            ea = ea,
            stratum = province,
            age = Q1_1,
            eversex = Q5_1A %>% na_if(0) %>% na_if(3) %>% `==`(1) %>% as.integer,
            sex12m = (eversex & (na_if(Q6_1, 3) == 1)) %>% as.integer,
            currpreg = if_else(eversex == 0 & (is.na(Q7_59)), FALSE,
                               Q7_59 == 1) %>% as.integer,
            hivstatus = factor(hivstat, 1:2, c("positive", "negative")),
            prev = as.integer(hivstatus == "positive"),
            weight = ibreal12,
            hivweight = ibreal1)


dat <- bind_rows(as_factor(za05),
                 as_factor(za08),
                 as_factor(za12))

dat$agegr <- cut(dat$age, 3:10*5, paste0(3:9*5, "-", 3:9*5+4), TRUE, FALSE)

des <- svydesign(~ea, strata=~surveyid+stratum,
                 data = dat %>% filter(!is.na(weight)) %>% mutate(hivstatus = "all"),
                 weights = ~weight,
                 nest=TRUE)

deshiv <- svydesign(~ea, strata=~surveyid+stratum,
                    data= dat %>% filter(!is.na(hivweight)),
                    weights = ~hivweight,
                    nest=TRUE)

svyprop <- function(formula, by, design){
  val <- left_join(svyby(formula, by, design, unwtd.count, na.rm=TRUE) %>%
                   rename(n = counts) %>% select(-se),
                   svyby(formula, by, design, svyciprop, vartype=c("se", "ci"), na.rm=TRUE))
  names(val) <- sub("^se\\..*", "se", names(val))
  val
}


za_preg <- bind_rows(
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr, des),
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr, deshiv),
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr,
          update(des, agegr = "15-49")),
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr,
          update(deshiv, agegr = "15-49"))
)

za_eversex <- bind_rows(
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr, des),
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr, deshiv),
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr,
          update(des, agegr = "15-49")),
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr,
          update(deshiv, agegr = "15-49"))
)

za_sex12m <- bind_rows(
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr, des),
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr, deshiv),
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr,
          update(des, agegr = "15-49")),
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr,
          update(deshiv, agegr = "15-49"))
)

za_pregprev <- bind_rows(
  svyprop(~prev, ~surveyid + country + survyear + agegr, subset(deshiv, currpreg == 1)),
  svyprop(~prev, ~surveyid + country + survyear + agegr,
          update(subset(deshiv, currpreg == 1), agegr = "15-49"))
)

za_prev <- bind_rows(
  svyprop(~prev, ~surveyid + country + survyear + agegr, deshiv),
  svyprop(~prev, ~surveyid + country + survyear + agegr, update(deshiv, agegr = "15-49"))
)


## Save results

save(za_preg, za_eversex, za_sex12m, za_pregprev, za_prev,
     file = "hsrc.rda")
