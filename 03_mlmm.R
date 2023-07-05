# clear workspace =====================================================================================================
rm(list = ls())

# set options =========================================================================================================
options(width = 120, digits = 8, scipen = 9, show.signif.stars = FALSE, mc.cores = 4)
runStan = TRUE # set to TRUE to run Stan models; if FALSE, will load saved results

buildTetrachoricPPMC = function(predictiveSamples, rawData, ...) {
  nObs = nrow(rawData)
  nItems = ncol(rawData)
  itemNames = names(rawData)
  
  # build posterior predictive samples matrix for tetrachoric correlations
  ttcPPMC = matrix(data = NA, nrow = nrow(predictiveSamples), ncol = nItems * (nItems - 1) / 2)
  tempNames = NULL
  for (item1 in 1:(nItems - 1)){
    for (item2 in (item1 + 1):nItems){
      tempNames = c(tempNames, paste0("TTC_", itemNames[item1], "_", itemNames[item2]))
    }
  }
  colnames(ttcPPMC) = tempNames
  simData = matrix(data = NA, nrow = nObs, ncol = nItems)
  
  pb = txtProgressBar()
  for (draw in 1:nrow(predictiveSamples)){
    for (item in 1:nItems){
      simData[, item] = predictiveSamples[draw, paste0("simY[", item, ",", 1:nObs, "]")]
    }
    
    col = 1
    for (item1 in 1:(nItems-1)){
      for (item2 in (item1+1):nItems){
        ttcPPMC[draw, col] = blatent:::findTetrachoric(data = simData, var1 = item1, var2 = item2, correct = .5)
        col = col + 1
      }
    }
    setTxtProgressBar(pb, draw/nrow(predictiveSamples))
  }
  
  # build data frame from draw
  col = 1
  obsTTC = matrix(data = NA, nrow = 1, ncol = nItems * (nItems - 1) / 2)
  for (item1 in 1:(nItems-1)){
    for (item2 in (item1+1):nItems){
      obsTTC[1, col] = blatent:::findTetrachoric(data = rawData, var1 = item1, var2 = item2, correct = .5)
      col = col + 1
    }
  }

  
  PPMCstats = NULL
  residTTC = matrix(data = 0, nrow = nItems, ncol = nItems)
  colnames(residTTC) = rownames(residTTC) = itemNames
  
  
  for (column in 1:length(obsTTC)){
    # build empirical distribution
    ttcdist = ecdf(ttcPPMC[, column])
    
    if (ttcdist(obsTTC[column]) > .95 || ttcdist(obsTTC[column]) < .05) {
      ttcFlag = TRUE
    } else {
      ttcFlag = FALSE
    }
    
    item1 = strsplit(x = colnames(ttcPPMC)[column], split = "_")[[1]][2]
    item2 = strsplit(x = colnames(ttcPPMC)[column], split = "_")[[1]][3]
    residTTCval = mean(ttcPPMC[, column]) - obsTTC[column]
    residTTC[item1, item2] = residTTC[item2, item1] = residTTCval
    
    # tabulate item correlation model fit statistics
    PPMCstats = rbind(PPMCstats, 
                      data.frame(
                        item1 = item1,
                        item2 = item2,
                        obsTTC = obsTTC[column],
                        estTTCmean = mean(ttcPPMC[, column]),
                        estTTCsd = sd(ttcPPMC[, column]),
                        estTTC05 = quantile(ttcPPMC[, column], probs = 0.05),
                        estTTC95 = quantile(ttcPPMC[, column], probs = 0.95),
                        residTTC = residTTCval,
                        absResidTTC = abs(mean(ttcPPMC[, column]) - obsTTC[column]),
                        PPPvalue = ttcdist(obsTTC[column]),
                        ttcFlag = ttcFlag
                      )
    )
    rownames(PPMCstats) = NULL
  }
  
  close(pb)
  
  # order by discrepancy 
  PPMCstats = PPMCstats[order(PPMCstats$absResidTTC, decreasing = TRUE), ]
  return(list(PPMCstats = PPMCstats, residTTC = residTTC))
}

# Package installation ================================================================================================
needed_packages = c("ggplot2", "cmdstanr", "HDInterval", "bayesplot", "loo")
for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage == FALSE) {
    install.packages(needed_packages[i])
  }
  library(needed_packages[i], character.only = TRUE)
}

# =====================================================================================================================
# Data Import, Manipulation, and Description
load("data/modelingData.RData")

correctReponseItems = names(modelingData)[grep(x = names(modelingData), pattern = "score")]
correctResponseData = modelingData[correctReponseItems]

nItems = length(grep(x = names(modelingData), pattern = "score"))
nPersons = nrow(modelingData)
nSchools = length(unique(modelingData$schoolID))
school = modelingData$schoolID


# =====================================================================================================================
# estimate multilevel multivariate two-level empty model

if (runStan){
  # compile stan model
  model2LE_Stan =  cmdstan_model(stan_file = "Stan/03_model01.stan")

  # build stan data file
  model2LE_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    priorItemRIsdMean = rep(0, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorRIcorrLKJparam = 1.0
  )

  # estimate model with MCMC
  model2LE_Samples = model2LE_Stan$sample(
    data = model2LE_Data,
    seed = 0606202301,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model2LE_Summary = model2LE_Samples$summary()
  model2LE_Draws = model2LE_Samples$draws()
  model2LE_Loo = model2LE_Samples$loo(variables = "personLike")  

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples2LE = model2LE_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc2LE = buildTetrachoricPPMC(predictiveSamples2LE, rawData)

  save(model2LE_Summary, model2LE_Draws, model2LE_Loo, ppmc2LE, file = "models/03_model01.RData")

} else {
  load("models/03_model01.RData")
}

# Assess convergence: summary of all parameters
max(
  model2LE_Summary$rhat[
    grep(pattern = "gamma|itemRIsd|itemRIcorr", x = model2LE_Summary$variable)
    ],
  na.rm = TRUE  
)


# parameter results
print(
  model2LE_Summary[
    grep(pattern = "gamma|itemRIsd|itemRIcorr\\[", x = model2LE_Summary$variable),
    ], 
  n = Inf
)

print(
  model2LE_Summary[
    grep(pattern = "itemRI\\[62", x = model2LE_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmc2LE$PPMCstats)
heatmap(ppmc2LE$residTTC)

# =====================================================================================================================
# adding person latent variable (without item discrimination/factor loading)

if (runStan){

  # compile stan model
  model2LWP_Stan =  cmdstan_model(stan_file = "Stan/03_model02.stan")

  # build stan data file
  model2LWP_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    priorItemRIsdMean = rep(0, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorRIcorrLKJparam = 1.0,
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1
  )

  model2LWP_Samples = model2LWP_Stan$sample(
    data = model2LWP_Data,
    seed = 0606202302,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000, 
    refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model2LWP_Summary = model2LWP_Samples$summary()
  model2LWP_Draws = model2LWP_Samples$draws()
  model2LWP_Loo = model2LWP_Samples$loo(variables = "personLike")

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples2LWP = model2LWP_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc2LWP = buildTetrachoricPPMC(predictiveSamples2LWP, rawData)

  save(model2LWP_Summary, model2LWP_Draws, model2LWP_Loo, ppmc2LWP, file = "models/03_model02.RData")

} else {
  load("models/03_model02.RData")
}

# Assess convergence: summary of all parameters
max(
  model2LWP_Summary$rhat[
    grep(pattern = "gamma|itemRI|personTheta", x = model2LWP_Summary$variable)
    ],
  na.rm = TRUE  
)


# parameter results
print(
  model2LWP_Summary[
    grep(pattern = "personTheta", x = model2LWP_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmc2LWP$PPMCstats)
heatmap(ppmc2LWP$residTTC)

 # compare relative model fit
loo_compare(
  list(
    model2LE = model2LE_Loo, 
    model2LWP = model2LWP_Loo
  )
)

# =====================================================================================================================
# adding person latent variable (with item discrimination/factor loading)

if (runStan){

  # compile stan model
  model2LWPD_Stan =  cmdstan_model(stan_file = "Stan/03_model03.stan")

  # build stan data file
  model2LWPD_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 1 * diag(nItems),
    school = school,
    priorItemRIsdMean = rep(0, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorRIcorrLKJparam = 1.0,
    priorMeanWithinLambda = rep(0, nItems),
    priorCovWithinLambda = 1 * diag(nItems)
  )


  # create starting values for theta that keep direction of theta meaningful
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = mean)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore

  model2LWPD_Samples = model2LWPD_Stan$sample(
    data = model2LWPD_Data,
    seed = 0606202303,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000, 
    refresh = 100,
    init = function() list(theta = startTheta, 
                           withinLambda = rep(5, nItems))
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model2LWPD_Summary = model2LWPD_Samples$summary()
  model2LWPD_Draws = model2LWPD_Samples$draws()
  model2LWPD_Loo = model2LWPD_Samples$loo(variables = "personLike")

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples2LWPD = model2LWPD_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc2LWPD = buildTetrachoricPPMC(predictiveSamples2LWPD, rawData)

  save(model2LWPD_Summary, model2LWPD_Draws, model2LWPD_Loo, ppmc2LWPD, file = "models/03_model03.RData")
} else {
  load("models/03_model03.RData")
}

# assess convergence: summary of all parameters
max(
  model2LWPD_Summary$rhat[
    grep(pattern = "gamma|itemRI|personTheta|withinLambda|thetaSD", x = model2LWPD_Summary$variable)
    ],
  na.rm = TRUE  
)

# parameter results
print(
  model2LWPD_Summary[
    grep(pattern = "withinLambda", x = model2LWPD_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmc2LWPD$PPMCstats)
heatmap(ppmc2LWPD$residTTC)

# compare relative model fit
loo_compare(
  list(
    model2LE = model2LE_Loo, 
    model2LWP = model2LWP_Loo,
    model2LWPD = model2LWPD_Loo
  )
) 

# compare results with 2PL model
load("models/01_model02.RData")

personID = unlist(
  lapply(
    X = model2PL_Summary$variable[grep(pattern = "theta", x = model2PL_Summary$variable)],
    FUN = function(x){
      locOpen = gregexpr(pattern = "\\[", text = x)[[1]][1]
      locClose = gregexpr(pattern = "\\]", text = x)[[1]][1]
      return(as.numeric(substr(x, locOpen + 1, locClose - 1)))
    }
  )
)

thetaComparison = data.frame(
  person = personID,
  theta2PL_EAP = model2PL_Summary$mean[grep(pattern = "theta", x = model2PL_Summary$variable)],
  theta2LWPD_EAP = model2LWPD_Summary$mean[grep(pattern = "personTheta", x = model2LWPD_Summary$variable)]
)


plot(x = thetaComparison$theta2PL_EAP, y = thetaComparison$theta2LWPD_EAP, 
     xlab = "2PL EAP", ylab = "2LWPD EAP", 
     xlim = c(-3, 3), ylim = c(-3, 3))




# =====================================================================================================================
# adding person latent variable (with item discrimination/factor loading) with marker-item identification 
# to estimate a latent variable variance 

if (runStan){

  # compile stan model
  model2LWPDV_Stan =  cmdstan_model(stan_file = "Stan/03_model04.stan")

  # build stan data file
  model2LWPDV_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    priorItemRIsdMean = rep(0, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorRIcorrLKJparam = 1.0,
    priorMeanWithinLambda = rep(0, nItems - 1),
    priorCovWithinLambda = 1 * diag(nItems - 1),
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1
  )

  # create starting values for theta that keep direction of theta meaningful
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = sum)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore


  model2LWPDV_Samples = model2LWPDV_Stan$sample(
    data = model2LWPDV_Data,
    seed = 0606202304,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000, 
    refresh = 100,
    init = function() list(theta = startTheta, 
                           withinLambda = rep(5, nItems - 1),
                           personThetaSD = 1)
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model2LWPDV_Summary = model2LWPDV_Samples$summary()
  model2LWPDV_Draws = model2LWPDV_Samples$draws()
  model2LWPDV_Loo = model2LWPDV_Samples$loo(variables = "personLike")

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples2LWPDV = model2LWPDV_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc2LWPDV = buildTetrachoricPPMC(predictiveSamples2LWPDV, rawData)

  save(model2LWPDV_Summary, model2LWPDV_Draws, model2LWPDV_Loo, ppmc2LWPDV, file = "models/03_model04.RData")
} else {
  load("models/03_model04.RData")
}

# assess convergence: summary of all parameters
max(
  model2LWPDV_Summary$rhat[
    grep(pattern = "gamma|itemRI|personTheta|withinLambda|personThetaSD", x = model2LWPDV_Summary$variable)
    ],
  na.rm = TRUE  
)

# parameter results
print(
  model2LWPDV_Summary[
    grep(pattern = "itemRIsd", x = model2LWPDV_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmc2LWPDV$PPMCstats)
heatmap(ppmc2LWPDV$residTTC)

# compare relative model fit
loo_compare(
  list(
    model2LE = model2LE_Loo, 
    model2LWP = model2LWP_Loo,
    model2LWPD = model2LWPD_Loo,
    model2LWPDV = model2LWPDV_Loo
  )
)

# =====================================================================================================================
# Adding school-level theta (betweenTheta) with separate discrimination parameters from withinTheta

if (runStan){
  # compile stan model
  modelWPBS_Stan = cmdstan_model(stan_file = "Stan/03_model05.stan")

  # build stan data file
  modelWPBS_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    priorItemRIsdMean = rep(0, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorMeanWithinLambda = rep(0, nItems - 1),
    priorCovWithinLambda = 1 * diag(nItems - 1),
    priorMeanBetweenLambda = rep(0, nItems - 1),
    priorCovBetweenLambda = 1 * diag(nItems - 1),
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1,
    priorSchoolThetaSDmean = 0,
    priorSchoolThetaSDsd = 1
  )

  modelWPBS_Samples = modelWPBS_Stan$sample(
    data = modelWPBS_Data,
    seed = 0606202305,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000, 
    refresh = 100,
    init = function() list(withinLambda = rep(5, nItems - 1),
                           betweenLambda = rep(5, nItems - 1),
                           personThetaSD = 1,
                           schoolThetaSD = 1)
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelWPBS_Summary = modelWPBS_Samples$summary()
  modelWPBS_Draws = modelWPBS_Samples$draws()
  modelWPBS_Loo = modelWPBS_Samples$loo(variables = "personLike")

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamplesWPBS = modelWPBS_Samples$draws(variables = "simY", format = "draws_matrix")
  
  ppmcWPBS = buildTetrachoricPPMC(predictiveSamplesWPBS, rawData)

  save(modelWPBS_Summary, modelWPBS_Draws, modelWPBS_Loo, ppmcWPBS, file = "models/03_model05.RData")
} else {
  load("models/03_model05.RData")
}

# assess convergence: summary of all parameters
max(
  modelWPBS_Summary$rhat[
    grep(pattern = "gamma|itemRI|personTheta|withinLambda|betweenLambda|personThetaSD|schoolThetaSD", 
         x = modelWPBS_Summary$variable)
    ],
  na.rm = TRUE  
)

# parameter results
print(
  modelWPBS_Summary[
    grep(pattern = "personTheta\\[", 
         x = modelWPBS_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmcWPBS$PPMCstats)
heatmap(ppmcWPBS$residTTC)

# compare relative model fit
loo_compare(
  list(
    model2LE = model2LE_Loo, 
    model2LWP = model2LWP_Loo,
    model2LWPD = model2LWPD_Loo,
    model2LWPDV = model2LWPDV_Loo,
    modelWPBS = modelWPBS_Loo
  )
)

# =====================================================================================================================
# constrain within-school and between-school loadings to be equal

if (runStan){

  # compile stan model
  modelMLMMED_Stan = cmdstan_model(stan_file = "Stan/03_model06.stan")

  # build stan data file
  modelMLMMED_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    priorItemRIsdMean = rep(0, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorMeanLambda = rep(0, nItems - 1),
    priorCovLambda = 1 * diag(nItems - 1),
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1,
    priorSchoolThetaSDmean = 0,
    priorSchoolThetaSDsd = 1
  )

  modelMLMMED_Samples = modelMLMMED_Stan$sample(
    data = modelMLMMED_Data,
    seed = 0606202306,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000, 
    refresh = 100,
    init = function() list(lambda = rep(5, nItems - 1),
                           personThetaSD = 1,
                           schoolThetaSD = 1)
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelMLMMED_Summary = modelMLMMED_Samples$summary()
  modelMLMMED_Draws = modelMLMMED_Samples$draws()
  modelMLMMED_Loo = modelMLMMED_Samples$loo(variables = "personLike")

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamplesMLMMED = modelMLMMED_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmcMLMMED = buildTetrachoricPPMC(predictiveSamplesMLMMED, rawData)

  save(modelMLMMED_Summary, modelMLMMED_Draws, modelMLMMED_Loo, ppmcMLMMED, file = "models/03_model06.RData")

} else {
  load("models/03_model06.RData")
}

# assess convergence: summary of all parameters
max(
  modelMLMMED_Summary$rhat[
    grep(pattern = "gamma|lambda|personTheta|personThetaSD|schoolThetaSD", 
         x = modelMLMMED_Summary$variable)
    ],
  na.rm = TRUE  
)

# parameter results
print(
  modelMLMMED_Summary[
    grep(pattern = "itemRIsd", 
         x = modelMLMMED_Summary$variable),
    ], 
  n = Inf
)



# parameter results
print(
  modelMLMMED_Summary[
    grep(pattern = "gamma|lambda|personTheta|personThetaSD|schoolThetaSD", 
         x = modelMLMMED_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmcMLMMED$PPMCstats)
heatmap(ppmcMLMMED$residTTC)

# compare relative model fit
loo_compare(
  list(
    model2LE = model2LE_Loo, 
    model2LWP = model2LWP_Loo,
    model2LWPD = model2LWPD_Loo,
    model2LWPDV = model2LWPDV_Loo,
    modelWPBS = modelWPBS_Loo,
    modelMLMMED = modelMLMMED_Loo
  )
)

# =====================================================================================================================
# removing item-level random intercepts (forcing homogeneous item intercepts across schools)

if (runStan){

  modelMLMMEDnoRI_Stan = cmdstan_model(stan_file = "Stan/03_model07.stan")

  # build stan data file
  modelMLMMEDnoRI_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    priorMeanLambda = rep(0, nItems - 1),
    priorCovLambda = 1 * diag(nItems - 1),
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1,
    priorSchoolThetaSDmean = 0,
    priorSchoolThetaSDsd = 1
  )

  modelMLMMEDnoRI_Samples = modelMLMMEDnoRI_Stan$sample(
    data = modelMLMMEDnoRI_Data,
    seed = 0606202307,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000, 
    refresh = 100,
    init = function() list(lambda = rep(5, nItems - 1),
                           personThetaSD = 1,
                           schoolThetaSD = 1)
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelMLMMEDnoRI_Summary = modelMLMMEDnoRI_Samples$summary()
  modelMLMMEDnoRI_Draws = modelMLMMEDnoRI_Samples$draws()
  modelMLMMEDnoRI_Loo = modelMLMMEDnoRI_Samples$loo(variables = "personLike")

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamplesMLMMEDnoRI = modelMLMMEDnoRI_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmcMLMMEDnoRI = buildTetrachoricPPMC(predictiveSamplesMLMMEDnoRI, rawData)

  save(
    modelMLMMEDnoRI_Summary, 
    modelMLMMEDnoRI_Draws,
    modelMLMMEDnoRI_Loo, 
    ppmcMLMMEDnoRI, 
    file = "models/03_model07.RData"
  )

} else {
  load("models/03_model07.RData")
}

# assess convergence: summary of all parameters
max(
  modelMLMMEDnoRI_Summary$rhat[
    grep(pattern = "gamma|lambda|personTheta|personThetaSD|schoolThetaSD", 
         x = modelMLMMEDnoRI_Summary$variable)
    ],
  na.rm = TRUE  
)

# parameter results
print(
  modelMLMMEDnoRI_Summary[
    grep(pattern = "gamma|lambda|personThetaSD|schoolThetaSD", 
         x = modelMLMMEDnoRI_Summary$variable),
    ], 
  n = Inf
)

# check model fit compared with 2PL model
# show model fit results
View(ppmcMLMMEDnoRI$PPMCstats)
heatmap(ppmcMLMMEDnoRI$residTTC)

# compare relative model fit
loo_compare(
  list(
    model2LE = model2LE_Loo, 
    model2LWP = model2LWP_Loo,
    model2LWPD = model2LWPD_Loo,
    model2LWPDV = model2LWPDV_Loo,
    modelWPBS = modelWPBS_Loo,
    modelMLMMED = modelMLMMED_Loo,
    modelMLMMEDnoRI = modelMLMMEDnoRI_Loo
  )
)

# =====================================================================================================================
# predict theta at both levels with frlunch

if (runStan){

  modelMLMMEDlunch_Stan =  cmdstan_model(stan_file = "Stan/03_model08.1.stan") 

  # build stan data file
  modelMLMMEDlunch_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 10 * diag(nItems),
    school = school,
    frlunch = modelingData$frlunch,
    priorItemRIsdMean = rep(-5, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorMeanLambda = rep(0, nItems - 1),
    priorCovLambda = 1 * diag(nItems - 1),
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1,
    priorSchoolThetaSDmean = 0,
    priorSchoolThetaSDsd = 1,
    priorFrlunchThetaSDmean = 0,
    priorFrlunchThetaSDsd = 1,
    priorFrlunchRIsdMean = 0,
    priorFrlunchRIsdSD = 1,
    priorSlopeSchoolThetaFRlunchMean = 0,
    priorSlopeSchoolThetaFRlunchSD = 1,
    priorSlopePersonThetaFRlunchMean = 0,
    priorSlopePersonThetaFRlunchSD = 1
  )

    # create starting values for theta that keep direction of theta meaningful
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = sum)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore

  modelMLMMEDlunch_Samples = modelMLMMEDlunch_Stan$sample(
    data = modelMLMMEDlunch_Data,
    seed = 0606202308,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 25,
    iter_sampling = 100, 
    refresh = 10,
    init = function() list(lambda = rep(5, nItems - 1),
                           personTheta = startTheta)
  )

  mcmc_trace(modelMLMMEDlunch_Samples$draws(variables = c("schoolRIsd", "schoolRIcov", "schoolRIcorr")))
  mcmc_trace(modelMLMMEDlunch_Samples$draws(variables = c("schoolThetaSD")))

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelMLMMEDlunch_Summary = modelMLMMEDlunch_Samples$summary()
  modelMLMMEDlunch_Draws = modelMLMMEDlunch_Samples$draws()
  modelMLMMEDlunch_Diag = modelMLMMEDlunch_Samples$sampler_diagnostics()

  save(
    modelMLMMEDlunch_Summary, 
    modelMLMMEDlunch_Draws,
    modelMLMMEDlunch_Diag,
    file = "models/03_model08.RData"
  )

} else {
  load("models/03_model08.RData")
}

# assess convergence: summary of all parameters
max(
  modelMLMMEDlunch_Summary$rhat[
    grep(pattern = "gamma|lambda|personTheta|personThetaSD|schoolThetaSD|frlunchThetaSD|frlunchRIsd|slopeSchoolThetaFRlunch|slopePersonThetaFRlunch", 
         x = modelMLMMEDlunch_Summary$variable)
    ],
  na.rm = TRUE  
)

View(
  modelMLMMEDlunch_Summary[
    grep(pattern = "slope|personThetaSD|schoolThetaSD|frlunchRIsd|lambda|gamma", 
         x = modelMLMMEDlunch_Summary$variable),
    ]
)

mcmc_trace(modelMLMMEDlunch_Samples$draws(), pars = paste0("gamma[", 1:9, "]"))

# parameter results
print(
  modelMLMMEDlunch_Summary[
    grep(pattern = "slope", 
         x = modelMLMMEDlunch_Summary$variable),
    ], 
  n = Inf
)



# =====================================================================================================================
# predict theta at person level with frlunch and allow school level theta to covary with frlunch random intercept

if (runStan){


  modelMLMMED2lunch_Stan =  cmdstan_model(stan_file = "Stan/03_model09.stan") 

  # build stan data file
  modelMLMMED2lunch_Data = list(
    nPersons = nPersons,
    nSchools = nSchools,
    nItems = nItems,
    Y = t(correctResponseData),
    priorMeanGamma = rep(0, nItems),
    priorCovGamma = 1 * diag(nItems),
    school = school,
    frlunch = modelingData$frlunch,
    priorItemRIsdMean = rep(-5, nItems),
    priorItemRIsdSD = rep(1, nItems),
    priorMeanLambda = rep(0, nItems - 1),
    priorCovLambda = 1 * diag(nItems - 1),
    priorPersonThetaSDmean = 0,
    priorPersonThetaSDsd = 1,
    priorSlopePersonThetaFRlunchMean = 0,
    priorSlopePersonThetaFRlunchSD = 1,
    priorSchoolRIsdMean = rep(0, 2),
    priorSchoolRIsdSD = rep(1, 2),
    priorSchoolRIcorrLKJparam= 1.0
  )

    # create starting values for theta that keep direction of theta meaningful
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = sum)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore

  modelMLMMED2lunch_Samples = modelMLMMED2lunch_Stan$sample(
    data = modelMLMMED2lunch_Data,
    seed = 0606202308,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 25,
    iter_sampling = 10000, 
    refresh = 10,
    save_warmup = TRUE,
    init = function() list(lambda = rep(5, nItems - 1),
                            personTheta = startTheta,
                            schoolRIcovL = diag(2),
                            schoolRIsd = rep(1, 2))
  )

# View(
#   modelMLMMED2lunch_Samples$summary(variables = c("schoolRIsd", "schoolRIcov", "schoolRIcorr"))
# )

# View(
#   modelMLMMED2lunch_Samples$summary(variables = c("lambda"))
# )

View(
  modelMLMMED2lunch_Samples$summary(variables = c("schoolPseudoR2"))
)

# mcmc_trace(modelMLMMED2lunch_Samples$draws(variables = c("schoolRIsd", "schoolRIcov", "schoolRIcorr")))
# mcmc_trace(modelMLMMED2lunch_Samples$draws(variables = c("gamma")))

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
modelMLMMED2lunch_Summary = modelMLMMED2lunch_Samples$summary()
modelMLMMED2lunch_Draws = modelMLMMED2lunch_Samples$draws()
modelMLMMED2lunch_Diag = modelMLMMED2lunch_Samples$sampler_diagnostics()

  save(
    modelMLMMED2lunch_Summary, 
    modelMLMMED2lunch_Draws,
    modelMLMMED2lunch_Diag,
    file = "models/03_model09.RData"
  )

} else {
  load("models/03_model09.RData")
}

# assess convergence: summary of all parameters
max(
  modelMLMMEDlunch_Summary$rhat[
    grep(pattern = "gamma|lambda|personTheta|personThetaSD|schoolThetaSD|frlunchThetaSD|frlunchRIsd|slopeSchoolThetaFRlunch|slopePersonThetaFRlunch", 
         x = modelMLMMEDlunch_Summary$variable)
    ],
  na.rm = TRUE  
)

mcmc_trace(modelMLMMEDlunch_Samples$draws(), pars = paste0("gamma[", 1:9, "]"))

# parameter results
print(
  modelMLMMED2lunch_Summary[
    grep(pattern = "gamma|lambda|personThetaSD|schoolThetaSD|frlunchThetaSD|frlunchRIsd|slopeSchoolThetaFRlunch|slopePersonThetaFRlunch", 
         x = modelMLMMEDlunch_Summary$variable),
    ], 
  n = Inf
)
