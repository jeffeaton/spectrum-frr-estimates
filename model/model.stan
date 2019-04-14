data {
  int<lower=0> N;
  int<lower=0> Kall;
  int<lower=0> Khivp;

  matrix[N, Kall] Xall;
  matrix[N, Khivp] Xhivp;

  vector[N] y;
  vector[N] pys;

  vector[N] hivstatus;
  matrix[N, 8] hivdist; // weighted HIV population on ART
  int n_country;
  int countryidx[N];
  int n_agegr;
  int agegridx[N];
  vector[6] log_frr_cd4_mean;
  matrix[6, 6] log_frr_cd4_cov;
}
parameters {
  vector[Kall] beta_all;
  vector[Khivp] beta_hivp;
  vector[n_agegr] lfrr_art;
  row_vector[6] log_frr_cd4;
  vector[n_country] u_country;
  real<lower = 0> sigma_u_country;
}
model {

  // informative prior for CD4 frr
  log_frr_cd4 ~ multi_normal(log_frr_cd4_mean, log_frr_cd4_cov);

  // diffuse prior for other coefficients
  beta_all ~ normal(0, 5);
  lfrr_art ~ normal(0, 5);
  beta_hivp ~ normal(0, 5);

  // country-level random effects for HIV FRR
  u_country ~ normal(0, sigma_u_country);
  sigma_u_country ~ normal(0, 1);


  // likelihood
  {
    vector[N] lfx;
    vector[N] lfall;
    vector[N] lfhivp;
    vector[N] fhivp;
    matrix[N, 8] frr_all;
    vector[N] lfrr;
    
    lfall = Xall * beta_all;
    
    if(Khivp > 0)
      lfhivp = Xhivp * beta_hivp;
    else
      lfhivp = rep_vector(0.0, N);
    
    if(n_country > 0)
    lfall = lfall + u_country[countryidx] .* hivstatus;
    
    frr_all[ , 1] = rep_vector(1.0, N) .* exp(lfhivp);
    frr_all[ , 2:7] = rep_matrix(exp(log_frr_cd4), N) .* rep_matrix(exp(lfhivp), 6);
    frr_all[ , 8] = exp(lfrr_art[agegridx]);
    
    lfrr = log(rows_dot_product(hivdist, frr_all));
    
    lfx = lfall + hivstatus .* lfrr;
    target += y .* lfx - exp(lfx) .* pys - log(pys);  // poisson with non-integer counts
  }
}
