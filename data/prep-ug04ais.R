library(tidyverse)
library(rdhs)
library(survey)

ugir <- read_zipdata("~/Documents/Data/DHS/Uganda/2004-05/UGIR5IDT.ZIP", readfn=haven::read_dta)
ugar <- read_zipdata("~/Documents/Data/DHS/Uganda/2004-05/UGAR5IDT.ZIP", readfn=haven::read_dta)

ug04 <- merge(ugir, ugar, by.x=c("v001", "v002", "v003"), by.y=c("hivclust", "hivnumb", "hivline"), all.x=TRUE)

ug04 <- ug04 %>%
  filter(v012 %in% 15:49) %>%
  mutate(surveyid = "UG2004AIS",
         survyear = 2004,
         country = "Uganda",
         ea = v001,
         stratum = paste(as_factor(v024), as_factor(v025)),
         age = v012,
         agegr = cut(age, 3:10*5, paste0(3:9*5, "-", 3:9*5+4), TRUE, FALSE),
         eversex = v525 == 0,
         sex12m = v766b > 0,
         currpreg = as.logical(v213),
         prev = as.integer(na_if(hiv03, 7) > 0),
         hivstatus = factor(prev, 0:1, c("negative", "positive")),
         weight = v005 / 1e6,
         hivweight = hiv05 / 1e6)


des <- svydesign(~ea, strata=~surveyid+stratum,
                 data = ug04 %>% filter(!is.na(weight)) %>% mutate(hivstatus = "all"),
                 weights = ~weight,
                 nest=TRUE)

deshiv <- svydesign(~ea, strata=~surveyid+stratum,
                    data = ug04 %>% filter(!is.na(hivweight)),
                    weights = ~hivweight,
                    nest=TRUE)

svyprop <- function(formula, by, design){
  val <- left_join(svyby(formula, by, design, unwtd.count, na.rm=TRUE) %>%
                   rename(n = counts) %>% select(-se),
                   svyby(formula, by, design, svyciprop, vartype=c("se", "ci"), na.rm=TRUE))
  names(val) <- sub("^se\\..*", "se", names(val))
  val
}


ug04_preg <- bind_rows(
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr, des),
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr, deshiv),
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr, update(des, agegr = "15-49")),
  svyprop(~currpreg, ~surveyid + country + survyear + hivstatus + agegr, update(deshiv, agegr = "15-49"))
)

ug04_eversex <- bind_rows(
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr, des),
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr, deshiv),
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr, update(des, agegr = "15-49")),
  svyprop(~eversex, ~surveyid + country + survyear + hivstatus + agegr, update(deshiv, agegr = "15-49"))
)

ug04_sex12m <- bind_rows(
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr, des),
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr, deshiv),
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr, update(des, agegr = "15-49")),
  svyprop(~sex12m, ~surveyid + country + survyear + hivstatus + agegr, update(deshiv, agegr = "15-49"))
)


ug04_pregprev <- bind_rows(
  svyprop(~prev, ~surveyid + country + survyear + agegr, subset(deshiv, currpreg == 1)),
  svyprop(~prev, ~surveyid + country + survyear + agegr,
          update(subset(deshiv, currpreg == 1), agegr = "15-49"))
)

ug04_prev <- bind_rows(
  svyprop(~prev, ~surveyid + country + survyear + agegr, deshiv),
  svyprop(~prev, ~surveyid + country + survyear + agegr, update(deshiv, agegr = "15-49"))
)


## Save results

save(ug04_preg, ug04_eversex, ug04_sex12m, ug04_pregprev, ug04_prev,
     file = here::here("data", "ug04ais.rda"))
