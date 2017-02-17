#!/usr/bin/env Rscript
library(rstan)
library(parallel)
options(mc.cores = parallel::detectCores())

load("dat/shw_scale.Rdata")

# set up an environmental data frame for stan
# start with just cells where we have all data
# include linear terms plus interactions for region
rows <- complete.cases(shw_scale)
env_dat <- with(shw_scale[rows,], data.frame(intercept=1,
                                      bathy=bathy,
                                      chla=chla,
                                      sst=sst,
                                      varsst=varsst,
                                      dist=dist,
                                      east_region = as.integer(region=='E')))
env_dat <- within(env_dat, {
  east_chla <- east_region*chla
  east_sst <- east_region*sst
  east_varsst <- east_region*varsst
  east_dist <- east_region*dist
})

# set up remaining stan data
stanDat <- list(
  num_cells = nrow(env_dat),
  num_env_vars = ncol(env_dat),
  env = as.matrix(env_dat),
  gps_counts = shw_scale[rows,'ninds'],
  # counts = shw_scale[rows,'count']
  abundance = shw_scale[rows,'count']
)

# launch stan
mod <- stan('shear.stan', data=stanDat, chains=4, iter=1000, cores=4)

