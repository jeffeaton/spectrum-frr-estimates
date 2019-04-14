library(tidyverse)
library(here)

load(here("data", "dhs.rda"))
load(here("data", "kais.rda"))
load(here("data", "hsrc.rda"))
load(here("data", "phia.rda"))
load(here("data", "mics.rda"))
load(here("data", "ug04ais.rda"))


tfr <- grep("\\_tfr", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows

asfr <- grep("\\_asfr", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows

currpreg <- grep("\\_preg$", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows

eversex <- grep("\\_eversex", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows

sex12m <- grep("\\_sex12m", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows

pregprev <- grep("\\_pregprev", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows %>%
  rename(pregprev = prev)

prev <- grep("\\_prev", ls(), value=TRUE) %>%
  lapply(get) %>%
  lapply(type.convert) %>%
  bind_rows 

lst <- list(tfr = tfr,
            asfr = asfr,
            currpreg = currpreg,
            eversex = eversex,
            sex12m  =  sex12m,
            pregprev = pregprev,
            prev = prev)

lst <- lst %>%
  lapply(function(x){if(!exists("tips", x)) x$tips <- 0; x}) %>%
  lapply(mutate,
         year = survyear - tips,
         type = sub(".{6}(.*)", "\\1", surveyid))

alldat <- lst %>%
  Map("[[<-", ., "indicator", value = names(.)) %>%
  Map(function(x, s){names(x) <- sub(s, "est", names(x)); x}, ., names(.)) %>%
  bind_rows

#' Assign region

reg <- c("Lesotho" = "South",
         "Mozambique" = "South",
         "Namibia" = "South",
         "South Africa" = "South",
         "Eswatini" = "South",
         "Zimbabwe" = "South",
         "Burundi" = "East",
         "Ethiopia" = "East",
         "Kenya" = "East",
         "Malawi" = "East",
         "Rwanda" = "East",
         "Tanzania" = "East",
         "Uganda" = "East",
         "Zambia" = "East",
         "Angola" = "West/Central",
         "Benin" = "West/Central",
         "Burkina Faso" = "West/Central",
         "Cameroon" = "West/Central",
         "Chad" = "West/Central",
         "Congo" = "West/Central",
         "Congo Democratic Republic" = "West/Central",
         "Cote d'Ivoire" = "West/Central",
         "Gabon" = "West/Central",
         "Gambia" = "West/Central",
         "Ghana" = "West/Central",
         "Guinea" = "West/Central",
         "Liberia" = "West/Central",
         "Mali" = "West/Central",
         "Niger" = "West/Central",
         "Nigeria" = "West/Central",
         "Senegal" = "West/Central",
         "Sierra Leone" = "West/Central",
         "Togo" = "West/Central")

alldat <- alldat %>%
  mutate(country = fct_recode(country, "Eswatini" = "Swaziland")) %>%
  right_join(data.frame(country = names(reg), region = reg), .) %>%
  filter(indicator != "eversex",
         indicator != "sex12m" | agegr == "15-19" & hivstatus == "all",
         indicator != "prev" | agegr == "15-49",
         indicator != "pregprev" | agegr == "15-49")
  

#' Create model dataset

moddat <- alldat %>%
  filter(indicator %in% c("currpreg", "asfr"),
         agegr != "15-49",
         hivstatus != "all",
         tips < 3,
         !is.na(region)) %>%
  rename(outcome = indicator) %>%
  mutate(n = if_else(outcome == "currpreg", as.numeric(n), pys),
         x = n * est,
         tips = if_else(outcome == "currpreg", -1, tips),
         year = survyear - tips) %>%
  select(surveyid, country, region, survyear, hivstatus, agegr, tips, year, type, outcome, est, se, x, n) %>%
  mutate(yearbefore = factor(tips + 1),
         hivstatus = as.integer(hivstatus == "positive"),
         currpreg = as.integer(outcome == "currpreg")) %>%
  left_join(
    eversex %>%
    filter(agegr == "15-19", hivstatus == "all") %>%
    mutate(eversex  = 100*eversex) %>%
    select(surveyid, eversex15to19 = eversex)
  ) %>%
  left_join(
    sex12m %>%
    filter(agegr == "15-19", hivstatus == "all") %>%
    mutate(sex12m  = 100*sex12m) %>%
    select(surveyid, sex12m)
  )


#' Distribution by CD4 stage and on ART from Spectrum 2018

pop1 <- readRDS(here::here("data", "pop1.rds"))

pop1w <- pop1 %>%
  filter(hivstatus == "positive") %>%
  mutate(cd4art = fct_relevel(cd4art, "500pl", "350to500", "250to350", "200to250", "100to200", "50to100", "bel50", "art6mospl")) %>%
  group_by(country, year, agegr, cd4art) %>%
  summarise(prop = sum(prop)) %>%
  spread(cd4art, prop, fill = 0) %>%
  as.data.frame

moddat <- moddat %>% left_join(pop1w)

save(alldat, moddat, pop1w, file = here::here("data", "data.rda"))
