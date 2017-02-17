data
{
  // sample sizes
  int <lower=0> nCounts; // counts ONLY
  int <lower=0> nCountsGPS; // counts AND GPS
  int <lower=0> nGPS; // GPS only
  int <lower=0> nEnvVars;
  
  // environmental matrices
  matrix [nCounts, nEnvVars] envCount;
  matrix [nCountsGPS, nEnvVars] envCountGPS;
  matrix [nGPS, nEnvVars] envGPS;
  
  // response variables
  int <lower=0> counts [nCounts]; // count data where we have NO gps data
  int <lower=0> countsGPS [nCountsGPS]; // count data where we also have gps data
  int <lower=0> gpsCounts [nCountsGPS]; // GPS counts where we also have aerial count data
  int <lower=0> gps [nGPS]; // GPS counts where we have no aerial survey data
  
  // offsets
  vector <lower=0> [nCounts] surveyEffortCounts;
  vector <lower=0> [nCountsGPS] surveyEffortCountsGPS;
}
transformed data
{
  matrix [nCounts + nCountsGPS + nGPS, nEnvVars] env; // all environmental variables smashed together
  
  env[1:nCounts,] = envCount;
  env[(nCounts+1):(nCounts+nCountsGPS), ] = envCountGPS;
  env[(nCounts + nCountsGPS + 1):(nCounts + nCountsGPS + nGPS),] = envGPS;
}
parameters
{
  vector [nEnvVars] beta;

  // latent expected GPS counts
  vector [nCountsGPS] gpsCountsRateLog;
  vector [nGPS] gpsRateLog;
  
  // because we cannot estimate the "effort" for GPS, we have an offset parameter that we
  // will learn from the data (capitalizing on the fact that the GPS and count data overlap)
  // we also need a variance hyperparameter, because we can't assume GPS tracks scale directly
  // to the count data
  real <lower=0> gpsScale;
  real <lower=0> gpsSD;
}
transformed parameters
{
  vector [nCounts + nCountsGPS + nGPS] abunRateLog; // log expected "true" abundance
  // vector [num_gps] abun_rate_gps_log;
  // vector [num_gps] gps_mean_log;
  // 
  abunRateLog = env * beta;
  // countRateLog = abunRateLog + log(countEffort)
  // abun_rate_gps_log = env_gps * beta;
  // gps_mean_log = abun_rate_gps_log + log(gps_scale);
}
model
{
  
  // count observation model; for now we assume our counts are direct observations of abundance
  // for the future, it might be a good idea to challenge this assumption
  // abundance ~ some_function(counts);
  
  // abundance model
  // abundance is a partially latent variable; it is observed via counts and (imperfectly) via GPS tracks
  counts ~ poisson_log(abunRateLog[1:nCounts] + log(surveyEffortCounts));
  countsGPS ~ poisson_log(abunRateLog[(nCounts+1):(nCounts+nCountsGPS)] + log(surveyEffortCountsGPS));

  // GPS MODEL
  // observation layer
  // gps counts are based centered around the mean abundance, plus a correction
  // to scale the gps counts to the true abundance
  // the expectation of the gps counts is the scale*abundance_rate
  // we model the log of this hierarchically with a gaussian prior
  gpsCounts ~ poisson_log(gpsCountsRateLog);
  gps ~ poisson_log(gpsRateLog);
  
  // hierarchical layer
  // the log rate for the GPS counts is itself a random variable centered around the scaled latent abundance
  gpsCountsRateLog ~ normal(abunRateLog[(nCounts+1):(nCounts+nCountsGPS)] + log(gpsScale), gpsSD);
  gpsRateLog ~ normal(abunRateLog[(nCounts+nCountsGPS+1):(nCounts+nCountsGPS+nGPS)] + log(gpsScale), gpsSD);
 
  
  // PRIORS for hyperparameters
  gpsSD ~ cauchy(0,10);
  gpsScale ~ cauchy(0,10);
  beta ~ cauchy(0,5);
  // 
}
