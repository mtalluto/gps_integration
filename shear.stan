data
{
  int <lower=0> num_cells;
  int <lower=0> num_env_vars;
  matrix environment [num_cells, num_env_vars];
  int gps_counts [num_cells];
  int counts [num_cells];
}
parameters
{
  vector <lower=0> abundance [num_cells];
  vector beta [num_env_vars];
  real <lower=0> gps_scale;
  real <lower=0> gps_sd;
}
model
{
  vector abun_rate [num_cells];
  vector gps_rate [num_cells];
  vector gps_mean [num_cells];
  
  // count observation model; for now we assume our counts are direct observations of abundance
  // for the future, it might be a good idea to challenge this assumption
  abundance = counts;
  
  // abundance model
  // abundance is a partially latent variable; it is observed via counts and (imperfectly) via GPS tracks
  abundance ~ poisson_log(abun_rate);
  abun_rate = beta * environment;
  
  // GPS MODEL
  // observation layer
  // gps counts are based centered around the mean abundance, plus a correction
  // to scale the gps counts to the true abundance
  // the expectation of the gps counts is the scale*abundance_rate
  // we model the log of this hierarchically with a gaussian prior
  gps_counts ~ poisson_log(gps_rate);
  gps_mean = abun_rate + log(gps_scale);
  gps_rate ~ normal(gps_mean, gps_sd);
  
  
  // PRIORS for hyperparameters
  gps_sd ~ cauchy(0,10);
  gps_scale ~ cauchy(0,10);
  beta ~ cauchy(0,5);
  
}