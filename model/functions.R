
log_frr_cd4_mean <- c(-0.103740703947107, -0.269890076205486, -0.582569410767653, 
                      -1.010608650165, -1.31446137695435, -1.42969881877763)

log_frr_cd4_cov <- matrix(c(0.0104778941807238, 0.00450151854182414, 0.00263954171411713, 
                            0.00154030864509876, 0.00224333672835132, 0.0036528378566598, 
                            0.00450151854182414, 0.0174458631891387, 0.00190050733332657, 
                            -0.00231858301282934, -0.00303277163476953, -0.00299485520149612, 
                            0.00263954171411713, 0.00190050733332657, 0.0358413960600581, 
                            0.00372620187623964, -0.00463894165922082, -0.00644900671576809, 
                            0.00154030864509876, -0.00231858301282934, 0.00372620187623964, 
                            0.0725647052225069, 0.0469542674869795, 0.0400435167052158, 0.00224333672835132, 
                            -0.00303277163476953, -0.00463894165922082, 0.0469542674869795, 
                            0.211131249777309, 0.206195413828095, 0.0036528378566598, -0.00299485520149612, 
                            -0.00644900671576809, 0.0400435167052158, 0.206195413828095, 
                            0.405721174766334),
                          nrow=6, ncol=6)

frr_cd4 <- c(1.0, 0.91, 0.77, 0.57, 0.38, 0.29, 0.28)


make_sdat <- function(f = ~-1 + currpreg + surveyid:agegr + yearbefore:hivstatus:agegrhiv,
                      fhivp = ~-1 + region:agegrhiv:hivstatus + age15to19hiv:sex12m_z,
                      dat){
  
  ## drop unused factor levels
  dat <- droplevels(dat)
  
  hivdist <- dat[c("500pl", "350to500", "250to350", "200to250", 
                   "100to200", "50to100", "bel50", "art6mospl")]
  
  sex12m_zmean <- 40
  sex12m_zsd <- 15
  dat$sex12m_z <- (dat$sex12m - sex12m_zmean) / sex12m_zsd

  mf <- model.frame(f, dat)
  mt <- attr(mf, "terms")
  Xall <- model.matrix(mt, mf)

  Xall <- Xall[ , !grepl(paste0("yearbefore", levels(dat$yearbefore)[1]), colnames(Xall))]

  mfhivp <- model.frame(fhivp, dat)
  mthivp <- attr(mfhivp, "terms")
  Xhivp <- model.matrix(mthivp, mfhivp)

  ## Check GLM fit to ensure no dropped levels
  if(ncol(Xhivp))
    glmfit <- summary(glm(x ~ -1 + Xall + Xhivp, "quasipoisson", dat, offset=log(n)))
  else
    glmfit <- summary(glm(x ~ -1 + Xall, "quasipoisson", dat, offset=log(n)))
  if(any(glmfit$aliased))
    stop("non-identified model formula")
  
  sdat <- list(N = nrow(Xall),
               Kall = ncol(Xall),
               Khivp = ncol(Xhivp),
               Xall = Xall,
               Xhivp = Xhivp,
               y = dat$x,
               pys = dat$n,
               hivdist = hivdist,
               hivstatus = dat$hivstatus,
               frr_cd4 = frr_cd4,
               log_frr_cd4_mean = log_frr_cd4_mean,
               log_frr_cd4_cov = log_frr_cd4_cov,
               n_country = length(levels(factor(dat$country))),
               countryidx = as.integer(factor(dat$country)),
               countryf = factor(dat$country),
               n_agegr = length(levels(dat$agegrhiv)),
               agegridx = as.integer(factor(dat$agegrhiv)),
               agegrf = factor(dat$agegrhiv),
               dat = dat,
               mf = mf,
               terms = mt,
               sex12m_zmean = sex12m_zmean,
               sex12m_zsd = sex12m_zsd)

  sdat
}
