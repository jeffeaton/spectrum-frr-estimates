library(tidyverse)

asfrcols <- c("surveyid", "country", "survyear", "hivstatus", "agegr", 
              "tips", "births", "pys", "asfr", "se")
tfrcols <- c("surveyid", "country", "survyear", "hivstatus", "tips", 
             "tfr", "se")
pregcols <- c("surveyid", "country", "survyear", "hivstatus", "agegr", 
              "n", "currpreg", "se", "ci_l", "ci_u")

phiayear <- c(Malawi = 2016, Zambia = 2016, Zimbabwe = 2016, Swaziland = 2017,
              Tanzania = 2017, Lesotho = 2017, Uganda = 2017)
phiaid <- setNames(paste0(c(Malawi = "MW", Zambia = "ZM", Zimbabwe = "ZW", Swaziland = "SZ",
                            Tanzania = "TZ", Lesotho = "LS", Uganda = "UG"),
                          phiayear, "PHIA"), names(phiayear))


phia_asfr <- read_csv("~/Documents/Data/PHIA/CDC/fertility/fertility-by-hiv.csv") %>%
  rename(asfr = fert_rate) %>%
  mutate(surveyid = phiaid[country],
         survyear = phiayear[country],
         hivstatus = c(all="all", hivn="negative", hivp="positive")[hivstatus],
         tips = 0,
         asfr = asfr / 1e3,
         ci_l = ci_l / 1e3,
         ci_u = ci_u / 1e3,
         se = (ci_u - ci_l) / (2*qnorm(0.975))) %>%
  filter(!is.na(surveyid)) %>%
  select(asfrcols)
         
phia_tfr <- phia_asfr %>%
  group_by(surveyid, country, survyear, hivstatus, tips) %>%
  summarise(tfr = 5*sum(asfr),
            se = 5*sqrt(sum(se^2)))


####  Currently pregnant by HIV status  ####

preg <- read_csv("~/Documents/Data/PHIA/CDC/fertility/pregnant-by-hiv.csv")

phia_preg <- preg %>%
  rename(currpreg = preg) %>%
  mutate(surveyid = phiaid[country],
         survyear = phiayear[country],
         agegr = sub("Total", "15-49", agegr),
         currpreg = currpreg / 100,
         ci_l = ci_l / 100,
         ci_u = ci_u / 100,
         se = (ci_u - ci_l) / (2 * qnorm(0.975))
         ) %>%
  filter(!is.na(surveyid)) %>%
  select(pregcols)
  


## Ever had sex

evsx <- read_csv("~/Documents/Data/PHIA/CDC/fertility/eversex-sex12m-modcontr.csv") %>%
  mutate(surveyid = phiaid[country],
         country = country,
         survyear = phiayear[country],
         agegr = sub("Total", "15-49", agegr),
         est = est / 100,
         ci_l = ci_l / 100,
         ci_u = ci_u / 100,
         se = (ci_u - ci_l) / (2 * qnorm(0.975))
         ) %>%
  filter(!is.na(surveyid)) %>%
  select(surveyid, country, survyear, hivstatus, agegr, n, est, se, ci_l, ci_u, indicator, sex)

phia_eversex <- evsx %>% filter(indicator == "eversex", sex == "female") %>%
  rename(eversex = est) %>%
  select(-indicator, -sex)

#' Note: sex12m is conditional on eversex. So pool exp(log(eversex) + log(sex12m))
phia_sex12m <- evsx %>%
  filter(indicator %in% c("eversex", "sex12m"), sex == "female") %>%
  group_by(surveyid, country, survyear, hivstatus, agegr) %>%
  summarise(n = max(n),
            sex12m = exp(sum(log(est))),
            se = sum((indicator == "eversex") * se),
            ci_l = sex12m - qnorm(0.975) * se,
            ci_u = sex12m + qnorm(0.975) * se) 

## Prevalence among pregnant women
## We are hacking this a bit...

phia_pregprev <- phia_preg %>% filter(agegr == "15-49") %>%
  select(-n, -se, -ci_l, -ci_u) %>%
  spread(hivstatus, currpreg) %>%
  mutate(prev = positive * (negative - all) / (all * (negative - positive))) %>%
  left_join(preg %>%
            filter(hivstatus == "all", agegr == "Total") %>%
            select(country, n = x)
            ) %>%
  mutate(se = sqrt(1.5 * prev * (1 - prev) / n),  # assume DEFF = 1.5
         ci_l = prev - qnorm(0.975) * se,
         ci_u = prev + qnorm(0.975) * se) %>%
  select(-all, -negative, -positive)

save(phia_preg, phia_tfr, phia_asfr, phia_eversex, phia_sex12m, phia_pregprev,
     file = here::here("data", "phia.rda"))
