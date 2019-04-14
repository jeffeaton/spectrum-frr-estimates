library(tidyverse)
library(rstan)
library(ggsci)
library(ggrepel)
library(gridExtra)


#' # Figure 1: summary
#'

load(here::here("data", "data.rda"))

table(alldat$country) %>% length

table(alldat$type)

surveys <- moddat %>%
  mutate(outcome = if_else(outcome == "asfr", paste0(outcome, tips), outcome)) %>%
  group_by(surveyid, country, survyear, region, type, outcome) %>%
  summarise(n = sum(n)) %>%
  ungroup %>%
  spread(outcome, n) %>%
  droplevels %>%
  mutate(country2 = sub("Congo Democratic Republic", "DR Congo", country),
         country2 = sub("Cote d'Ivoire", "Côte d'Ivoire", country2),
         region2 = fct_recode(region, "Southern" = "South")) %>%
  arrange(region2, country2) %>%
  mutate(region = fct_inorder(region),
         region2 = fct_inorder(region2),
         country = fct_inorder(country),
         country2 = fct_inorder(country2),
         type = fct_relevel(type, "DHS", "AIS", "PHIA", "KAIS", "HSRC"))


fig1data <- surveys %>%
  mutate(countryidx = as.integer(fct_rev(country2)) + c(0, 1, 2)[match(region2, c("West/Central", "Southern", "East"))],
         data = if_else(is.na(asfr0), "currpreg only", "birth history"))

country_labels <- fig1data %>%
  select(region, region2, country, country2, countryidx) %>% unique

fig1data %>%
  ggplot(aes(survyear, countryidx, color = type, shape = data, size = currpreg)) +
  annotate("text", 2001.5, c(9, 21.5, 29.5), label = c("West/Central", "Southern", "Eastern"), angle = 90, fontface = "bold") +
  geom_hline(yintercept = c(18, 25)) +
  geom_point(stroke = 1) +
  scale_x_continuous("survey year", breaks = 2003:2017, minor_breaks = NULL) +
  scale_y_continuous(element_blank(), breaks = country_labels$countryidx, minor_breaks = NULL,
                     labels = country_labels$country2,
                     position = "right", expand = expand_scale(add = 0.6)) +
  scale_size_area("sample size", c(3000, 9000, 15000), max_size = 3) +
  scale_color_nejm() +
    scale_shape_manual(values = c("birth history" = 19, "currpreg only" = 1)) + 
  theme_light(10) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         size = guide_legend(order = 3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing = unit(0, "cm"),
        legend.key.size = unit(0, "pt"),
        legend.justification = c(1, 0),
        legend.text.align=0,
        plot.margin = margin(t = 1, l = 15, unit = "pt")) +
  coord_cartesian(xlim = c(2003, 2017), clip = "off")

ggsave("figure1.png", w = 82, h = 150, units = "mm")


#' Analysis for text

surveys %>% count(country)
surveys %>% count(type)
surveys %>% count(is.na(asfr0))
surveys %>% count(is.na(currpreg))

surveys %>% count(is.na(asfr1), type)

surveys %>% filter(is.na(asfr2), asfr0)
surveys %>% filter(type == "AIS")

surveys %>% filter(survyear >= 2015)
surveys %>% filter(survyear < 2015)

table(surveys$type)

#' First paragraph of results

surveys %>% count(wt = currpreg)

moddat %>% count(outcome, wt = n)
moddat %>% count(outcome, hivstatus, wt = n)

moddat %>% count(outcome, wt = x)
moddat %>% count(outcome, hivstatus, wt = x)


#' # Table 1: estimates for model parameters

f1 <- readRDS(here::here("model", "fits", "fit_all.rds"))
f2 <- readRDS(here::here("model", "fits", "fit_bef15.rds"))
f3 <- readRDS(here::here("model", "fits", "fit_fert.rds"))
f4 <- readRDS(here::here("model", "fits", "fit_currpreg.rds"))

fits <- list("1: all" = f1, "2: bef2015" = f2, "3: fert only" = f3, "4: currpreg only" = f4)

tab1dat <- fits %>%
  lapply(function(x) data.frame(param = c(colnames(x$dat$Xhivp),
                                          paste0("log_frr_cd4_", c("350-499", "250-349", "200-249", "100-199", "50-99", "<50")),
                                          paste0("lfrr_art_", c("15-19", "20-24", "25-29", "30-34", "35-49"))),
                                summary(x$fit, c("beta_hivp", "log_frr_cd4", "lfrr_art"))$summary,
                                check.names = FALSE)) %>%
  Map("[[<-", ., "model", value = names(.)) %>%
  Map("[[<-", ., "modelid", value = 1:4) %>%
  bind_rows()

tab1 <- tab1dat %>%
  mutate(value = sprintf("%.3f (%.3f)", mean, sd)) %>%
  select(param, model, value) %>%
  spread(model, value)

tab1 %>% write.csv("table1.csv", row.names = FALSE)

round(f1$dat$log_frr_cd4_mean, 2) %>% t %>% t

sprintf("%.2f (%.2f)", f1$dat$log_frr_cd4_mean, sqrt(diag(f1$dat$log_frr_cd4_cov))) %>% t %>% t

#' For text

tab1dat %>%
  mutate(agegrhiv = sub(".*agegrhiv([^:]+):.*", "\\1", param)) %>%
  filter(agegrhiv == "15-19", modelid == 1) %>%
  transmute(param,
            est = exp(mean) %>% round(2),
            ci_l = exp(`2.5%`) %>% round(2),
            ci_u = exp(`97.5%`) %>% round(2))

tab1dat %>%
  filter(param == "age15to19hiv:sex12m_z", modelid == 1) %>%
  transmute(param,
            est = exp(mean) %>% round(2), 
            ci_l = exp(`2.5%`) %>% round(2),
            ci_u = exp(`97.5%`) %>% round(2))

tab1dat %>%
  mutate(agegrhiv = sub(".*agegrhiv([^:]+):.*", "\\1", param)) %>%
  filter(agegrhiv %in% c("20-24", "25-29"), modelid == 1) %>%
  transmute(param,
            est = exp(mean) %>% round(2),
            ci_l = exp(`2.5%`) %>% round(2),
            ci_u = exp(`97.5%`) %>% round(2))

tab1dat %>%
  mutate(agegrhiv = sub(".*agegrhiv([^:]+):.*", "\\1", param)) %>%
  filter(agegrhiv %in% c("30-34", "35-49"), modelid == 1) %>%
  transmute(param,
            est = exp(mean) %>% round(2),
            ci_l = exp(`2.5%`) %>% round(2),
            ci_u = exp(`97.5%`) %>% round(2))

tab1dat %>%
  filter(grepl("^lfrr\\_art", param), modelid == 1) %>%
  transmute(param,
            est = exp(mean) %>% round(2),
            ci_l = exp(`2.5%`) %>% round(2),
            ci_u = exp(`97.5%`) %>% round(2))


#' # Table 2: estimates for country-specific effects

tab2dat <- fits %>%
  lapply(function(x) data.frame(country = c(levels(x$dat$countryf), "sigma_u"),
                                summary(x$fit, c("u_country", "sigma_u_country"))$summary,
                                check.names = FALSE)) %>%
  Map("[[<-", ., "model", value = names(.)) %>%
  Map("[[<-", ., "modelid", value = 1:4) %>%
  bind_rows()

tab2dat %>%
  mutate(value = sprintf("%.3f (%.3f)", mean, sd)) %>%
  select(country, model, value) %>%
  spread(model, value) %>%
  right_join(country_labels %>% select(region2, country2, country), .) %>%
  arrange(region2, country2) %>%
  write.csv("table2.csv", na = "", row.names = FALSE)

tab2dat %>%
  filter(modelid == 1) %>%
  transmute(country,
            est = exp(mean) %>% round(2),
            ci_l = exp(`2.5%`) %>% round(2),
            ci_u = exp(`97.5%`) %>% round(2))


#' # Figure 2: Summary of estimates

leg_theme <- theme_light(10) +
  theme(legend.text = element_text(size = 8),
                   legend.key.size = unit(12, "pt"),
                   legend.title = element_text(size = 8, face = "bold"),
                   plot.tag = element_text(size = 12, face = "bold"))

## A: Age effects

fig2adat <- f1 %>%
  {
    data.frame(param = colnames(.$dat$Xhivp),
               t(as.matrix(.$fit, "beta_hivp")),
               check.names = FALSE)
  } %>%
  mutate(region = sub("region([^:]+):.*", "\\1", param),
         agegrhiv = sub(".*:agegrhiv([^:]+):.*", "\\1", param)) %>%
  gather(sampleid, value, -param, -region, -agegrhiv)

fig2a <- fig2adat %>%
  filter(param != "age15to19hiv:sex12m_z") %>%
  mutate(region = fct_recode(region, "West/\nCentral" = "West/Central")) %>%
  ggplot(aes(agegrhiv, value, fill = region)) +
  geom_hline(yintercept = 0.0, linetype="dashed") +
  geom_violin(linetype="blank",
              position=position_dodge(width=0.8), alpha = 0.9) +
  stat_summary(fun.y=mean, geom="point", size=1.5, color=1,
               position=position_dodge(width=0.8), show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous("FRR untreated / ART <6mos\nCD4 >500 (log scale)",
                     breaks = -2:2 * log(1.5),
                     labels = function(x) round(exp(x), 2)) +
  xlab(element_blank()) +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  guides(fill = guide_legend(ncol = 2, direction = "horizontal",
                             legend.key.height = unit(12, "pt"),
                             legend.text = element_text(size = 3))) +
  labs(tag = "A") +
  leg_theme +
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.box.background = element_rect(color = "grey"),
        legend.key.size = unit(13, "pt"),
        ## legend.title = element_text(size = 8),
        legend.text = element_text(size = 6.5))
  

#' B: CD4 effects

fig2bdat <- as.matrix(f1$fit, paste("log_frr_cd4")) %>%
  t %>%
  rbind(0, .) %>%
  data.frame(cd4 = c("≥500", "350-499", "250-349", "200-249", "100-199", "50-99", "<50"), ., check.names = FALSE) %>%
  mutate(cd4 = fct_inorder(cd4)) %>%
  gather(sampleid, value, -cd4)
               

fig2b <- fig2bdat %>%
  ggplot(aes(cd4, exp(value))) +
  geom_hline(yintercept = 1.0, linetype="dashed") +
  geom_violin(linetype="blank", fill="red3", alpha = 0.9) +
  stat_summary(fun.y=mean, geom="point", size=2, color=1,
               position=position_dodge(width=0.8)) +
  scale_y_continuous("Relative fertility", limits=c(0, 1.2)) +  
  xlab("CD4 category") +
  labs(tag = "B") +
  leg_theme +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))



#' C: ART effects

fig2cdat <- as.matrix(f1$fit, paste("lfrr_art")) %>%
  t %>%
  data.frame(agegrhiv = c("15-19", "20-24", "25-29", "30-34", "35-49"), check.names = FALSE) %>%
  mutate(agegrhiv = fct_inorder(agegrhiv)) %>%
  gather(sampleid, value, -agegrhiv)
               

fig2c <- fig2cdat %>%
  ggplot(aes(agegrhiv, value, fill = "ART >6mos")) +
  geom_hline(yintercept = 0.0, linetype="dashed") +
  geom_violin(linetype="blank", alpha = 0.9) +
  stat_summary(fun.y=mean, geom="point", size=2, color=1,
               position=position_dodge(width=0.8), show.legend = FALSE) +
  theme_light(9) +
  scale_y_continuous("\nFRR ART >6 mos (log scale)",
                     breaks = -2:2 * log(1.5),
                     labels = function(x) round(exp(x), 2)) +
  scale_fill_manual(element_blank(), values = c("ART >6mos" = "purple3")) +
  xlab(element_blank()) +
  labs(tag = "C") +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  leg_theme +
  theme(legend.position = c(0.99, 0.985),
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.box.background = element_rect(color = "grey"))


## Proportion sex12m among 15-19

nm <- colnames(f1$dat$Xhivp)
idx <- c(grep("^region.*:agegrhiv15-19:hivstatus$", nm), grep("sex12m", nm))
beta <- as.matrix(f1$fit, paste("beta_hivp"))[ , idx]
colnames(beta) <- nm[idx]

dat <- expand.grid(region = c("East", "South", "West/Centrall"),
                   agegrhiv = "15-19",
                   sex12m = seq(10, 70, 1)) %>%
  mutate(sex12m_z = (sex12m - f1$dat$sex12m_zmean) / f1$dat$sex12m_zsd)

fig2ddat <- data.frame(dat, model.matrix(~ -1 + region + sex12m_z, dat) %*% t(beta), check.names = FALSE) %>%
  gather(sampleid, value, -region, -agegrhiv, -sex12m, -sex12m_z) %>%
  mutate(evalue = exp(value))
  
fig2d <- fig2ddat %>%
  group_by(region, sex12m) %>%
  summarise(mean = mean(evalue),
            lower = quantile(evalue, 0.025),
            upper = quantile(evalue, 0.975)) %>%
  mutate(sex12m = sex12m / 100,
         x0 = if_else(sex12m == 0.40, sex12m, NA_real_)) %>%
  ggplot(aes(sex12m, mean, ymin = lower, ymax = upper, color = region, fill = region)) +
  geom_hline(yintercept = 1.0, linetype="dashed") +
  geom_ribbon(alpha=0.2, linetype = "blank") +
  geom_line(size=1) +
  geom_point(aes(x = x0), size = 2.0, show.legend = FALSE) +
  ylab("FRR age 15-19\nuntreated / ART <6mos") +
  scale_x_continuous("Proportion sexually active in past 12 months", labels = scales::percent_format(1)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(tag = "D") +
  leg_theme + 
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.box.background = element_rect(color = "grey"),
        legend.key.size = unit(13, "pt"),
        ## legend.title = element_text(size = 8),
        legend.text = element_text(size = 6.5))


fig2 <- arrangeGrob(fig2a, fig2b, fig2c, fig2d, ncol = 2)
ggsave("figure2.png", fig2, width = 173, height = 130, units = "mm")
  


#' # Figure 3: Country random effect

fig3data <- f1 %>%
  {data.frame(country = levels(.$dat$countryf),
              t(as.matrix(.$fit, "u_country")),
              check.names = FALSE)
  } %>%
  gather(sampleid, value, -country) %>%
  right_join(country_labels %>% select(region = region2, country2, country), .) %>%
  arrange(region, country2)

fig3 <- fig3data %>%
  ggplot(aes(country2, value, fill=region)) +
  geom_hline(yintercept = 0.0, linetype="dashed") +
  geom_violin(linetype="blank") +
  stat_summary(fun.y=mean, geom="point", size=2, color=1, show.legend = FALSE) +
  theme_light(10) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(element_blank()) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  ylab(expression(u[country])) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(12, "pt"),
        legend.title = element_text(size = 9, face = "bold"))

ggsave("figure3.png", fig3, width = 173, height = 70, units = "mm")


#' # Figure 4:

specbirths <- readRDS(here::here("data", "specbirths.rds"))

pop1 <- readRDS(here::here("data", "pop1.rds")) %>%
  filter(hivstatus == "positive") %>%
  mutate(cd4art = cd4art %>% fct_recode("≥500" = "500pl",
                                        "350-499" = "350to500",
                                        "250-349" = "250to350",
                                        "200-249" = "200to250",
                                        "100-199" = "100to200",
                                        "50-99" = "50to100",
                                        "<50" = "bel50") %>%
         fct_relevel("≥500", "350-499", "250-349", "200-249", "100-199", "50-99", "<50")) %>%
  group_by(country, year, agegr, cd4art) %>%
  summarise(prop = sum(prop))

system.time(
  agefrr_pred <- pop1 %>%
  left_join(
    alldat %>% filter(agegr == "15-19", hivstatus == "all", indicator == "sex12m") %>%
    group_by(country) %>%
    arrange(-survyear) %>%
    filter(row_number() == 1) %>%
    select(region, country, sex12m = est)
  ) %>%
  mutate(sex12m_z = (100 * sex12m - f1$dat$sex12m_zmean) / f1$dat$sex12m_zsd, 
         agegrhiv = fct_collapse(agegr, "35-49" = c("35-39", "40-44", "45-49")),
         age15to19hiv = as.integer(agegr == "15-19")) %>%
    left_join(fig3data %>% select(-region) %>% rename(u_country = value) %>% filter(sampleid %in% 1:2000)) %>%
    left_join(fig2adat %>% filter(param == "age15to19hiv:sex12m_z") %>% select(sampleid, beta_sex12m_z = value)) %>%
    left_join(fig2adat %>% filter(param != "age15to19hiv:sex12m_z") %>% rename(age_frr = value) %>% select(-param)) %>%
    left_join(fig2bdat %>% rename(cd4art = cd4, cd4_frr = value)) %>%
    left_join(fig2cdat %>% rename(art_frr = value)) %>%
    mutate(lfrr = u_country + if_else(cd4art == "art6mospl", art_frr, age_frr + cd4_frr + sex12m_z * age15to19hiv * beta_sex12m_z)) %>%
    group_by(region, country, year, agegr, sampleid) %>% 
    summarise(frr = sum(prop * exp(lfrr)))
)

asfrhiv_pred <- agefrr_pred %>%
  filter(country != "Botswana") %>%
  left_join(specbirths %>% filter(agegr != "15-49", outcome %in% c("births", "fpop", "hivprev")) %>% spread(outcome, value)) %>%
  mutate(asfr = births / fpop,
         births_hivp = births * frr * hivprev / (1 - (1 - frr) * hivprev),
         births_hivn = births - births_hivp,
         asfr_hivp = births_hivp / (fpop * hivprev),
         asfr_hivn = births_hivn / (fpop * (1 - hivprev)),
         preg_hivp = births_hivp * exp(-0.61705072),
         preg_hivn = births_hivn * exp(-0.61705072))


tfrhiv_pred <- asfrhiv_pred %>%
  group_by(region, country, year, sampleid) %>%
  summarise(tfr = 5 * sum(asfr),
            tfr_hivp = 5 * sum(asfr_hivp),
            tfr_hivn = 5 * sum(asfr_hivn),
            currpreg = sum(preg_hivp + preg_hivn) / sum(fpop),
            currpreg_hivp = sum(preg_hivp) / sum(fpop * hivprev),
            currpreg_hivn = sum(preg_hivn) / sum(fpop * (1.0 - hivprev)),
            hivprev = sum(hivprev * fpop) / sum(fpop),
            pregprev = sum(births_hivp) / sum(births)) %>%
  ungroup

fig4dat <- tfrhiv_pred %>%
  mutate(prev_ratio = pregprev / hivprev,
         tfr_ratio = tfr_hivp / tfr_hivn) %>%
  gather(key, value, -(region:sampleid)) %>%
  group_by(region, country, key, year) %>%
  summarise(mean = mean(value),
            median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  ungroup

figS1dat <- tfrhiv_pred %>%
  gather(key, value, -(region:sampleid)) %>%
  group_by(region, country, key, year) %>%
  summarise(mean = mean(value),
            median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  ungroup
            
saveRDS(tfrhiv_pred, "tfrhiv_pred.rds")
saveRDS(figS1dat, "figS1dat.rds")
saveRDS(fig4dat, "fig4dat.rds")

fig4dat <- readRDS("fig4dat.rds")

labc <- c("Eswatini", "Ethiopia", "Kenya", "Malawi", "Ghana", "DR Congo", "South Africa", "Senegal")


fig4a <- fig4dat %>%
  filter(key == "tfr_ratio",
         year %in% 1990:2020) %>%
  left_join(country_labels) %>%
  mutate(label = if_else(year == 2020 & country2 %in% labc, 
                         as.character(country2), NA_character_)) %>%
  ggplot(aes(year, mean, group = country, color = region2, label = label)) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  geom_line() +
  geom_text(nudge_x = 0.25, hjust = 0, size = 2.0, show.legend = FALSE) +
  scale_x_continuous(element_blank(), limits = c(1990, 2025)) +
  scale_y_continuous("TFR ratio:\nHIV positive / HIV negative women", limits = c(0.5, 1.25)) +
  scale_color_brewer("region", palette = "Set1") +
  labs(tag = "A") +
  leg_theme + 
  theme(legend.position = "none",
        legend.justification = c(1, 1),
        legend.box.background = element_rect(color = "grey"),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 6.5))

fig4b <- fig4dat %>%
  filter(key == "prev_ratio",
         year %in% 1990:2020) %>%
  left_join(country_labels) %>%
  mutate(label = if_else(year == 2020 & country2 %in% labc, 
                         as.character(country2), NA_character_)) %>%
  ggplot(aes(year, mean, group = country, color = region2, label = label)) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  geom_line() +
  geom_text(nudge_x = 0.25, hjust = 0, size = 2.0, show.legend = FALSE) +
  scale_x_continuous(element_blank(), limits = c(1990, 2025)) +
  scale_y_continuous("HIV prevalence ratio:\ncurrently pregnant / women 15-49", limits = c(0.35, 1.5)) +
  scale_color_brewer("region", palette = "Set1") +
  labs(tag = "B") +
  leg_theme + 
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.box.background = element_rect(color = "grey"),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 6.5))


fig4 <- arrangeGrob(fig4a, fig4b, ncol = 2)
ggsave("figure4.png", fig4, width = 173, height = 70, units = "mm")


#' Figure S1-S3

figS1dat <- readRDS("figS1dat.rds")

figS1dat <- figS1dat %>%
  mutate(hivstatus = factor(sub(".*[\\_](.*)", "\\1", key), c("hivp", "hivn"), c("positive", "negative")),
         key = sub("(.*)[\\_](.*)", "\\1", key)) %>%
  inner_join(country_labels)

plotdat <- alldat %>% filter(country %in% figS1dat$country,
                             (hivstatus != "all" | indicator %in% c("prev", "pregprev"))) %>%
  inner_join(country_labels)

figS1 <- figS1dat %>%
  filter(key == "tfr", !is.na(hivstatus), year %in% 2000:2018) %>%
  ggplot(aes(year, mean, ymin = lower, ymax = upper, color = hivstatus, fill = hivstatus)) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line() +
  geom_point(aes(year, est, color = hivstatus),
             plotdat %>% filter(indicator == "tfr", tips == 0),
             inherit.aes = FALSE,
             position = position_dodge(0.9),
             show.legend = FALSE) +
  geom_linerange(aes(year, est, ymin = est - qnorm(0.975)*se, ymax = est + qnorm(0.975)*se, color = hivstatus),
                 plotdat %>% filter(indicator == "tfr", tips == 0),
                 inherit.aes = FALSE,
                 position = position_dodge(0.9),
                 show.legend = FALSE) +
  facet_wrap(~country2) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(element_blank()) +
  scale_y_continuous("total fertiity rate (TFR)") +
  leg_theme +
  theme_light(8) +
  theme(legend.position = c(0.98, 0),
        legend.justification = c(1, 0),
        legend.direction = "horizontal")


figS2 <- figS1dat %>%
  filter(key == "currpreg", !is.na(hivstatus), year %in% 2000:2018) %>%
  ggplot(aes(year, mean, ymin = lower, ymax = upper, color = hivstatus, fill = hivstatus)) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line() +
  geom_point(aes(year, est, color = hivstatus),
             plotdat %>% filter(indicator == "currpreg", agegr == "15-49"),
             inherit.aes = FALSE,
             position = position_dodge(0.8),
             show.legend = FALSE) +
  geom_linerange(aes(year, est, ymin = ci_l, ymax = ci_u, color = hivstatus),
                 plotdat %>% filter(indicator == "currpreg", agegr == "15-49"),
                 inherit.aes = FALSE,
                 position = position_dodge(0.8),
                 show.legend = FALSE) +
  facet_wrap(~country2) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(element_blank()) +
  scale_y_continuous("Percent currently pregnant", 
                     labels = scales::percent_format(1)) +
  coord_cartesian(ylim = c(0.0, 0.18)) +
  leg_theme +
  theme_light(8) +
  theme(legend.position = c(0.98, 0),
        legend.justification = c(1, 0),
        legend.direction = "horizontal")


figS3surv <- plotdat %>%
  filter(indicator %in% c("pregprev", "prev"), agegr == "15-49") %>%
  mutate(population = factor(indicator, c("prev", "pregprev"), c("women 15-49", "pregnant women")))

figS3 <- figS1dat %>%
  filter(key %in% c("hivprev", "pregprev"), year %in% 2000:2018) %>%
  mutate(population = factor(key, c("hivprev", "pregprev"), c("women 15-49", "pregnant women"))) %>%
  ggplot(aes(year, mean, ymin = lower, ymax = upper, color = population, fill = population)) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line() +
   geom_point(aes(year, est, color = population),
              figS3surv,
              inherit.aes = FALSE,
              position = position_dodge(0.8),
              show.legend = FALSE) +
  geom_linerange(aes(year, est, ymin = ci_l, ymax = ci_u, color = population),
                 figS3surv,
                 inherit.aes = FALSE,
                 position = position_dodge(0.8),
                 show.legend = FALSE) +
  facet_wrap(~country2, scale = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(element_blank()) +
  scale_y_continuous("HIV prevalence", 
                     labels = scales::percent_format(0.1)) +
  leg_theme +
  theme_light(8) +
  theme(legend.position = c(0.98, 0),
        legend.justification = c(1, 0),
        legend.direction = "horizontal")

ggsave("figureS1.png", figS1, height = 7, width = 6, units = "in")
ggsave("figureS2.png", figS2, height = 7, width = 6, units = "in")
ggsave("figureS3.png", figS3, height = 6, width = 9.0, units = "in")


#' Summary parameter estimates

alldat %>% filter(agegr == "15-19", hivstatus == "all", indicator == "sex12m") %>%
    group_by(country) %>%
    arrange(-survyear) %>%
    filter(row_number() == 1) %>%
  select(region, country, survyear, surveyid, type, sex12m = est) %>%
    write.csv("sex12m.csv", row.names = FALSE)

fig3data %>%
  group_by(region, country) %>%
  summarise(est = mean(exp(value))) %>%
  write.csv("frr_country.csv", row.names = FALSE)

fig2adat %>%
  group_by(region, agegrhiv) %>%
  summarise(est = mean(exp(value)))  %>%
  write.csv("frr_age.csv", row.names = FALSE)

fig2bdat %>%
  group_by(cd4) %>%
  summarise(est = mean(exp(value))) %>%
  write.csv("frr_cd4.csv", row.names = FALSE)

fig2cdat %>%
  group_by(agegrhiv) %>%
  summarise(est = mean(exp(value))) %>%
  write.csv("frr_art.csv", row.names = FALSE)
