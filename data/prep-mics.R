library(tidyverse)
library(haven)
library(demogsurv)
library(survey)
library(rdhs)

## SB1       "Age at first sexual intercourse"                             
## SB2U      "Time since last sexual intercourse (unit)"                   
## SB2N      "Time since last sexual intercourse (number)"                 

## HA16               Tested for AIDS virus as part of antenatal care
## HA17              Received results from test during antenatal care
## HA18     Received consultation after testing during antenatal care
## HA20                         Tested for AIDS virus during delivery
## HA21                    Received results from test during delivery
## HA22             Tested for AIDS virus since test during pregnancy
## HA23                    Most recent time of testing for AIDS virus
## HA24                               Ever been tested for AIDS virus
## HA25                    Most recent time of testing for AIDS virus
## HA26                           Received results of AIDS virus test
## HA27                           Know a place to get AIDS virus test

## For MICS, ever tested in composite of HA16/HA17, HA20/HA21, and HA24/HA26.

set_na <- function(x, na_codes = 9){x[x %in% na_codes] <- NA; x}

svyprop <- function(formula, by, design){
  val <- left_join(svyby(formula, by, design, unwtd.count, na.rm=TRUE) %>%
                   rename(n = counts) %>% select(-se),
                   svyby(formula, by, design, svyciprop, vartype=c("se", "ci"), na.rm=TRUE))
  names(val) <- sub("^se\\..*", "se", names(val))
  val
}

do_mics <- function(zfile, country, cc){

  wm <- read_zipdata(zfile, "wm.sav", haven::read_spss) %>% as.data.frame()
  bh <- read_zipdata(zfile, "bh.sav", haven::read_spss) %>% as.data.frame()

  names(wm) <- tolower(names(wm))
  names(bh) <- tolower(names(bh))

  wm$survyear <- median(wm$wm6y)
  wm$surveyid <- paste0(cc, wm$survyear, "MICS")
  wm$country <- country

  print(paste(paste0(cc, median(wm$wm6y), "MICS"), median(wm$wm6y), country))

  ## Recode line number for MICS4
  if(is.null(bh$ln) & !is.null(bh$wm4))
    bh$ln <- bh$wm4

  nbh <- nrow(bh)

  wm$caseid <- seq_len(nrow(wm))
  bh <- merge(wm[c("caseid", "hh1", "hh2", "ln")], bh, all.y=TRUE)

  if(nrow(bh) > nbh)
    stop("you invented new births by merging wm and bh")

  ## Recode age for MICS4
  if(is.null(wm$wb2) & !is.null(wm$wm9))
    wm$wb2 <- wm$wm9
  
  ## Recode DOI variable for MICS4
  if(is.null(wm$wdoi) & !is.null(wm$cmcdoiw))
    wm$wdoi <- wm$cmcdoiw

  ## Recode chdob variable for MICS4
  if(is.null(bh$bh4c) & !is.null(bh$ccdob))
    bh$bh4c <- bh$ccdob
  
  
  wm$currpreg <- as.integer(set_na(wm$cp1, 9) == 1)

  if(!is.null(wm$sb1))
    wm <- wm %>% mutate(eversex = as.integer(set_na(sb1, c(97, 99)) > 0))

  if(!is.null(wm$sb2u))
    wm <- wm %>% mutate(sex12m = as.integer(eversex & (set_na(sb2u, 9) %in% 1:3)))
  else if(!is.null(wm$sb3u))
    wm <- wm %>% mutate(sex12m = as.integer(eversex & (set_na(sb3u, 9) %in% 1:3)))

  tfr <- calc_tfr(subset(wm, wm7 == 1), by=~surveyid+country+survyear,
                  tips = 0:6, bhdata = bh,
                  clusters = ~wm1, strata = NULL, dob = "wdob", intv = "wdoi",
                  weight = "wmweight", bvars = "bh4c")

  asfr <- calc_asfr(subset(wm, wm7 == 1), by=~surveyid+country+survyear,
                    tips = 0:6, bhdata = bh,
                    clusters = ~wm1, strata = NULL, dob = "wdob", intv = "wdoi",
                    weight = "wmweight", bvars = "bh4c", counts=TRUE)

  agegr <- 3:10*5
  wm$agegr <- cut(wm$wb2, agegr, demogsurv:::.epis_labels(agegr), TRUE, FALSE)
  des <- svydesign(~wm1, data=subset(wm, wb2 %in% 15:49), weights=~wmweight)

  if(!is.null(wm$eversex))
     eversex <- bind_rows(svyprop(~eversex, ~surveyid+country+survyear+agegr, des),
                          svyprop(~eversex, ~surveyid+country+survyear+agegr, update(des, agegr="15-49")))
  else
    eversex <- NULL

  if(!is.null(wm$sex12m))
    sex12m <- bind_rows(svyprop(~sex12m, ~surveyid+country+survyear+agegr, des),
                        svyprop(~sex12m, ~surveyid+country+survyear+agegr, update(des, agegr="15-49")))
  else
    sex12m <- NULL

  preg <- bind_rows(svyprop(~currpreg, ~surveyid+country+survyear+agegr, des),
                    svyprop(~currpreg, ~surveyid+country+survyear+agegr, update(des, agegr="15-49")))

  rownames(preg) <- rownames(eversex) <- rownames(tfr) <- rownames(asfr) <- NULL
  names(asfr) <- sub("^se[\\.\\_].*", "se", names(asfr))
  names(tfr) <- sub("^se[\\.\\_].*", "se", names(tfr))

  return(list(asfr = asfr, tfr = tfr, preg = preg, eversex = eversex, sex12m = sex12m))
}


########################
####  Analyse MICS  ####
########################

micsdir <- "~/Documents/Data/MICS/"
cbind(list.files(micsdir, "ZIP$", ignore.case=TRUE))

zfile <- file.path(micsdir, c("Malawi 2006 MICS_Datasets.zip",
                              "Malawi 2013-14 MICS_Datasets.zip",
                              "Swaziland 2010 MICS_Datasets.zip",
                              "Swaziland_MICS5_Datasets.zip",
                              "Zimbabwe 2009 MICS_Datasets.zip",
                              "Zimbabwe 2014 MICS_Datasets.zip"))
country <- c("Malawi", "Malawi", "Swaziland", "Swaziland", "Zimbabwe", "Zimbabwe")
cc <- c("MW", "MW", "SZ", "SZ", "ZW", "ZW")

mics <- Map(do_mics, zfile, country, cc)

mics <- do.call(Map, c(bind_rows, mics)) %>%
  lapply("[[<-", "hivstatus", value = "all")

Map(assign, paste0("mics_", names(mics)), mics, list(globalenv()))

mics_tfr <- mics_tfr %>% rename(se = se_tfr)
mics_asfr <- mics_asfr %>% rename(se = se_asfr)

save(mics_preg, mics_eversex, mics_sex12m, mics_tfr, mics_asfr,
     file = here::here("data", "mics.rda"))
