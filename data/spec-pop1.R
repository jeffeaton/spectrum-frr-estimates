library(tidyverse)
library(eppasm)

#########################################
####  Process Spectrum pop1 outputs  ####
#########################################

## Calculate distribution of women across CD4 and ART stages by age
## for reproductive aged women.


dir <- "~/Documents/Data/Spectrum files/2018 final/SSA"

pop1files <- bind_rows(
  data.frame(country = "Angola", file=list.files(dir, "Angola.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Botswana", file=list.files(dir, "Botswana.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Burkina Faso", file=list.files(dir, "Burkina.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Burundi", file=list.files(dir, "Burundi.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Cameroon", file=list.files(dir, "Cameroon.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Chad", file=list.files(dir, "Chad.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Congo", file=list.files(dir, "^Congo\\_.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Congo Democratic Republic", file=list.files(dir, "DRC.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Cote d'Ivoire", file=list.files(dir, "Ivoire.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Eswatini", file=list.files(dir, "Swaziland.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Ethiopia", file=list.files(dir, "Ethiopia.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Gabon", file=list.files(dir, "Gabon.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Gambia", file=list.files(dir, "Gambia.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Ghana", file=list.files(dir, "Ghana.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Guinea", file=list.files(dir, "^Guinee.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Kenya", file=list.files(dir, "Kenya.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Lesotho", file=list.files(dir, "Lesotho.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Liberia", file=list.files(dir, "Liberia.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Malawi", file=list.files(dir, "Malawi.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Mali", file=list.files(dir, "Mali.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Mozambique", file=list.files(dir, "MZ.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Niger", file=list.files(dir, "^Niger_.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Namibia", file=list.files(dir, "Namibia.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Rwanda", file=list.files(dir, "Rwanda.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Senegal", file=list.files(dir, "Senegal.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Sierra Leone", file=list.files(dir, "Sierra.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "South Africa", file=list.files(dir, "South Africa.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Tanzania", file=list.files(dir, "TZ.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Togo", file=list.files(dir, "Togo.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Uganda", file=list.files(dir, "Uganda.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Zambia", file=list.files(dir, "Zambia.*pop1\\.xlsx", full.names=TRUE)),
  data.frame(country = "Zimbabwe", file=list.files(dir, "^zim.*pop1\\.xlsx", full.names=TRUE))
)


read_pop1 <- function(pop1file, country, years = 1990:2021){

  print(country)

  pop1 <- lapply(as.character(years), readxl::read_excel, path = pop1file) %>%
    Map("[[<-", ., "year", value = years) %>%
    bind_rows %>%
    rename(sex = "..1",
           cd4 = "..2",
           artdur = `Single Age`) %>%
    filter(!is.na(sex)) %>%
    gather(age, pop, `0`:`80+`) %>%
    mutate(sex = fct_recode(sex, "male" = "Males", "female" = "Females"),
           age = sub("\\+", "", age) %>% type.convert,
           cd4 = sub("CD4 Count: ", "", cd4) %>% type.convert,
           artdur = sub("Duration: ", "", artdur) %>% type.convert,
           country = country) %>%
    filter(sex == "female",
           age %in% 15:49,
           cd4 != 0, artdur != 0,  # cd4 = 0 and artdur = 0 corresponde to total population estimates
           pop > 0)

  pop1 %>% as.data.frame
}           

pop1 <- parallel::mcMap(read_pop1, pop1file=pop1files$file, country=pop1files$country, mc.cores = parallel::detectCores())

pop1 <- bind_rows(pop1) %>%
  mutate(agegr = cut(age, 3:10*5, paste0(3:9*5, "-", 3:9*5+4), include.lowest=TRUE, right=FALSE),
         cd4 = factor(cd4, 1:8, c("hivn", "500pl", "350to500", "250to350", "200to250", "100to200", "50to100", "bel50")),
         artdur = factor(artdur, c(1:2, 6:8), c("hivn", "hivp", "art0mos", "art6mos", "art1yr")),
         hivstatus = if_else(cd4 == "hivn", "negative", "positive"),
         cd4art = if_else(artdur == "hivn", "hivn",
                  if_else(artdur %in% c("hivp", "art0mos"), as.character(cd4),
                  if_else(artdur %in% c("art6mos", "art1yr"), "art6mospl", NA_character_)))) %>%
  group_by(country, year, sex, agegr, hivstatus, cd4, artdur, cd4art) %>%
  summarise(pop = sum(pop)) %>%
  group_by(country, year, sex, agegr, hivstatus) %>%
  mutate(prop = pop / sum(pop)) %>%
  ungroup

saveRDS(pop1, here::here("data", "pop1.rds"))
        


#############################
####  Get births by age  ####
#############################

pjnz <- sub("\\_pop1.xlsx", "\\.PJNZ", basename(pop1files$file))
pop1files$pjnz <- sapply(pjnz, list.files, path=dir, full.names=TRUE, ignore.case = TRUE)

read_specbirths <- function(pjnz, country){

  print(country)

  demp <- read_specdp_demog_param(pjnz)
  sr <- read_hivproj_output(pjnz)

  agegrf <- cut(15:49, 3:10*5, paste0(3:9*5, "-", 3:9*5+4), TRUE, FALSE)
  mm <- model.matrix(~-1+agegr, data.frame(agegr=agegrf))
  dimnames(mm)[[2]] <- levels(agegrf)
  names(dimnames(mm))[2] <- "agegr"

  fpop <- sr$totpop[ 16:50, "Female", ]
  fpop <- (fpop[,-1] + fpop[,-ncol(fpop)]) / 2
  names(dimnames(fpop)) <- c("age", "year")
  births <- t(mm) %*% (demp$asfr[ , -1] * fpop)
  fpop <- t(mm) %*% fpop

  hivpop <- sr$hivpop[ 16:50, "Female", ]
  hivpop <- (hivpop[,-1] + hivpop[,-ncol(hivpop)]) / 2
  names(dimnames(hivpop)) <- c("age", "year")
  hivpop <- t(mm) %*% hivpop

  val <- list(births = reshape2::melt(births),
              fpop = reshape2::melt(fpop),
              hivpop = reshape2::melt(hivpop)) %>%
    Map("[[<-", ., "outcome", value = names(.)) %>%
    bind_rows %>%
    mutate(agegr = fct_expand(agegr, "15-49")) %>%
    bind_rows(
      data.frame(outcome = "births",
                 agegr = factor("15-49", levels(.$agegr)),
                 year = names(sr$births) %>% type.convert,
                 value = sr$births),
      data.frame(outcome = "hivpreg",
                 agegr = factor("15-49", levels(.$agegr)),
                 year = names(sr$hivpregwomen) %>% type.convert,
                 value = sr$hivpregwomen)
    )

  data.frame(country, val)
}


specbirths <- parallel::mcMap(read_specbirths, pop1files$pjnz, pop1files$country, mc.cores = parallel::detectCores())

specbirths <- specbirths %>%
  bind_rows %>%
  group_by(country, agegr, year, outcome) %>%
  summarise(value = sum(value)) %>%
  spread(outcome, value) %>%
  mutate(asfr = births / fpop,
         hivprev = hivpop / fpop) %>%
  gather(outcome, value, -(country:year)) %>%
  ungroup %>%
  filter(!is.na(value)) %>%
  as.data.frame

spectfr <- specbirths %>%
  filter(outcome == "asfr", agegr != "15-49") %>%
  group_by(country, year) %>%
  summarise(tfr = 5 * sum(value)) %>%
  filter(year %in% 2000:2020) %>%
  spread(year, tfr) %>%
  as.data.frame

saveRDS(specbirths, here::here("data", "specbirths.rds"))
saveRDS(spectfr, here::here("data", "spectfr.rds"))
