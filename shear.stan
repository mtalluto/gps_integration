data
{
  int <lower=0> num_cells;
  int <lower=0> num_env_vars;
  matrix [num_cells, num_env_vars] env;
  int gps_counts [num_cells];
  int abundance [num_cells];
  // int counts [num_cells];
}
parameters
{
  // we are not treating this as latent just yet
  // vector <lower=0> [num_cells] abundance;
  vector [num_env_vars] beta;
  vector [num_cells] gps_rate_log;
  real <lower=0> gps_scale;
  real <lower=0> gps_sd;
}
transformed parameters
{
  vector [num_cells] abun_rate_log;
  vector [num_cells] gps_mean_log;
  abun_rate_log = env * beta;

  gps_mean_log = abun_rate_log + log(gps_scale);
}
model
{
  // 
  // for(i in 1:num_cells) print(abun_rate_log[i]);
  
  // count observation model; for now we assume our counts are direct observations of abundance
  // for the future, it might be a good idea to challenge this assumption
  // abundance ~ some_function(counts);
  
  // abundance model
  // abundance is a partially latent variable; it is observed via counts and (imperfectly) via GPS tracks
  abundance ~ poisson_log(abun_rate_log);
  
  // GPS MODEL
  // observation layer
  // gps counts are based centered around the mean abundance, plus a correction
  // to scale the gps counts to the true abundance
  // the expectation of the gps counts is the scale*abundance_rate
  // we model the log of this hierarchically with a gaussian prior
  gps_counts ~ poisson_log(gps_rate_log);
  gps_rate_log ~ normal(gps_mean_log, gps_sd);
  
  
  // PRIORS for hyperparameters
  gps_sd ~ cauchy(0,10);
  gps_scale ~ cauchy(0,10);
  beta ~ cauchy(0,5);
  
}
