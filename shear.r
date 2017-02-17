#!/usr/bin/env Rscript
library("rstan")
library("parallel")
library("bayesplot")
library("raster")
options(mc.cores = parallel::detectCores())

load("dat/shw_scale.Rdata")

#### DELETE THIS LINE ONCE EFFORT IS ADDED
shw_scale$effort <- 1

# set up an environmental data frame for stan
# start with just cells where we have all data
# include linear terms plus interactions for region
envDatAll <- with(shw_scale, data.frame(intercept=1, bathy=bathy, chla=chla, sst=sst, varsst=varsst,
                  dist=dist, east_region = as.integer(region=='E')))

envDatAll <- within(envDatAll, {
  east_chla <- east_region*chla
  east_sst <- east_region*sst
  east_varsst <- east_region*varsst
  east_dist <- east_region*dist
})


# counts, no gps, and no environmental NAs
# counts, GPS, and no environmental NAs
# GPS only, no environmental NAs
countRows <- which(!is.na(shw_scale$count) & is.na(shw_scale$ninds) & complete.cases(envDatAll))
gpsCountRows <- which(!is.na(shw_scale$ninds) & !is.na(shw_scale$count) & complete.cases(envDatAll))
gpsRows <- which(!is.na(shw_scale$ninds) & is.na(shw_scale$count) & complete.cases(envDatAll))


# data for fitting
envCount <- envDatAll[countRows,]
envCountGPS <- envDatAll[gpsCountRows,]
envGPS <- envDatAll[gpsRows,]

# set up remaining stan data
stanDat <- list(
  nCounts = nrow(envCount),
  nCountsGPS = nrow(envCountGPS),
  nGPS = nrow(envGPS),
  nEnvVars = ncol(envCount),
  envCount = as.matrix(envCount),
  envCountGPS = as.matrix(envCountGPS),
  envGPS = as.matrix(envGPS),
  counts = shw_scale[countRows,'count'],
  countsGPS = shw_scale[gpsCountRows,'count'],
  gpsCounts = shw_scale[gpsCountRows,'ninds'],
  gps = shw_scale[gpsRows,'ninds'],
  surveyEffortCounts = shw_scale[countRows,'effort'],
  surveyEffortCountsGPS = shw_scale[gpsCountRows,'effort'])

envAbunOnly <- rbind(envCount, envCountGPS)
stanDatAbun <- list(nCounts = nrow(envAbunOnly),
   nEnvVars = ncol(envAbunOnly),
   envCount = as.matrix(envAbunOnly),
   counts = shw_scale[c(countRows, gpsCountRows),'count'],
   surveyEffortCounts = shw_scale[c(countRows, gpsCountRows),'effort'])

envGPSOnly <- rbind(envCountGPS, envGPS)
stanDatGPS <- list(nGPS = nrow(envGPSOnly),
  nEnvVars = ncol(envGPSOnly),
  envGPS = as.matrix(envGPSOnly),
  gps = shw_scale[c(gpsCountRows, gpsRows),'ninds'])

  
# launch stan
modFull <- stan('shear.stan', data=stanDat, chains=4, iter=1000, cores=4)
modAerial <- stan('shear_abun_only.stan', data=stanDatAbun, chains=4, iter=1000, cores=4)
modGPS <- stan('shear_gps_only.stan', data=stanDatGPS, chains=4, iter=1000, cores=4)

# visualize results
paramsFull <- as.array(modFull)
betasFull <- rstan::extract(modFull)$beta
paramsAerial <- as.array(modAerial)
betasAerial <- rstan::extract(modAerial)$beta
paramsGPS <- as.array(modGPS)
betasGPS <- rstan::extract(modGPS)$beta

# look at betas and hyperparams
# mcmc_dens(paramsFull, pars=c('gpsScale', 'gpsSD'), regex_pars='beta')

# predict abundance for all replicates across entire region
logAbunFull <- as.matrix(envDatAll) %*% t(betasFull)
logAbunAerial <- as.matrix(envDatAll) %*% t(betasAerial)
logAbunGPS <- as.matrix(envDatAll) %*% t(betasGPS)

# make posterior mean and SD rasters
abunFullMeanRas <- rasterFromXYZ(cbind(shw_scale$lon, shw_scale$lat, exp(rowMeans(logAbunFull))))
abunFullSDRas <- rasterFromXYZ(cbind(shw_scale$lon, shw_scale$lat, apply(logAbunFull, 1, sd)))
abunAerialMeanRas <- rasterFromXYZ(cbind(shw_scale$lon, shw_scale$lat, exp(rowMeans(logAbunAerial))))
abunAerialSDRas <- rasterFromXYZ(cbind(shw_scale$lon, shw_scale$lat, apply(logAbunAerial, 1, sd)))
abunGPSMeanRas <- rasterFromXYZ(cbind(shw_scale$lon, shw_scale$lat, exp(rowMeans(logAbunGPS))))
abunGPSSDRas <- rasterFromXYZ(cbind(shw_scale$lon, shw_scale$lat, apply(logAbunGPS, 1, sd)))

aCols <- colorRampPalette(rev(c('#032838', '#4E788F', '#80A4B7', '#BAD0D9', '#F1EBDC')), bias=3)(100)
pdf("results.pdf", h=10, w=10)
par(mfrow=c(3,2), mar=c(0.5,0.5,4,4))
plot(abunFullMeanRas, col=aCols, xaxt='n', yaxt='n', zlim=c(0, maxValue(abunFullMeanRas)), main="Mean Predicted Abundance (Integrated)")
plot(abunFullSDRas, col=heat.colors(100), xaxt='n', yaxt='n', zlim=c(0, maxValue(abunFullSDRas)), main="SD Predicted Abundance")
plot(abunAerialMeanRas, col=aCols, xaxt='n', yaxt='n', zlim=c(0, maxValue(abunAerialMeanRas)), main="Mean Predicted Abundance (Aerial Only)")
plot(abunAerialSDRas, col=heat.colors(100), xaxt='n', yaxt='n', zlim=c(0, maxValue(abunAerialSDRas)), main="SD Predicted Abundance")
plot(abunGPSMeanRas, col=aCols, xaxt='n', yaxt='n', zlim=c(0, maxValue(abunGPSMeanRas)), main="Mean Predicted Abundance (GPS Only)")
plot(abunGPSSDRas, col=heat.colors(100), xaxt='n', yaxt='n', zlim=c(0, maxValue(abunGPSSDRas)), main="SD Predicted Abundance")
dev.off()