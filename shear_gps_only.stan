data
{
  // sample sizes
  int <lower=0> nGPS;
  int <lower=0> nEnvVars;
  
  // environmental matrices
  matrix [nGPS, nEnvVars] envGPS;
  
  // response variables
  int <lower=0> gps [nGPS];
}
parameters
{
  vector [nEnvVars] beta;
}
transformed parameters
{
  vector [nGPS] abunRateLog; // log expected "true" abundance
  abunRateLog = envGPS * beta;
}
model
{
  // GPS MODEL
  gps ~ poisson_log(abunRateLog);
  
  beta ~ cauchy(0,5);
  // 
}
