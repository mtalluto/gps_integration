data
{
  int <lower=0> nCounts; // counts ONLY
  int <lower=0> nEnvVars;
  matrix [nCounts, nEnvVars] envCount;
  int counts [nCounts]; // count data where we have NO gps data
  vector <lower=0> [nCounts] surveyEffortCounts;
}
parameters
{
  vector [nEnvVars] beta;
}
transformed parameters
{
  vector [nCounts] abunRateLog;

  abunRateLog = envCount * beta;
}
model
{
  counts ~ poisson_log(abunRateLog + log(surveyEffortCounts));

  // PRIORS for hyperparameters
  beta ~ cauchy(0,5);
  
}
