# clear workspace =====================================================================================================
rm(list = ls())

# set options =========================================================================================================
options(width = 120, digits = 8, scipen = 9, show.signif.stars = FALSE, mc.cores = 4)

runStan = TRUE # set to TRUE to run Stan models; if FALSE, will load saved results

# install functions needed for analysis ===============================================================================
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
needed_packages = 
  c("ggplot2", "cmdstanr", "HDInterval", "bayesplot", "loo", "blatent", "utils", "mirt", "latticeExtra")
for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage == FALSE) {
    install.packages(needed_packages[i])
  }
  library(needed_packages[i], character.only = TRUE)
}

# load modeling data ==================================================================================================
load("data/modelingData.RData")


# determine data specs ================================================================================================
correctReponseItems = names(modelingData)[grep(x = names(modelingData), pattern = "score")]
correctResponseData = modelingData[correctReponseItems]
nItems = length(correctReponseItems)
nObs = nrow(correctResponseData)

# non-Bayesian analyses ===============================================================================================

# 1PL analysis with mirt
model1PL_MIRT = mirt(data = correctResponseData, model = 1, itemtype = "Rasch")
coef(model1PL_MIRT, IRTpars = FALSE)
coef(model1PL_MIRT, IRTpars = TRUE)
plot(model1PL_MIRT, type = "trace") 
plot(model1PL_MIRT, type = "infoSE") 
M2(model1PL_MIRT)

# 2PL analysis with mirt
model2PL_MIRT = mirt(data = correctResponseData, model = 1, itemtype = "2PL")
coef(model2PL_MIRT, IRTpars = FALSE)
coef(model2PL_MIRT, IRTpars = TRUE)
plot(model2PL_MIRT, type = "trace") 
plot(model2PL_MIRT, type = "infoSE") 
M2(model2PL_MIRT)

# model comparison
anova(model1PL_MIRT, model2PL_MIRT)

# Bayesian analyses ===================================================================================================

# 1PL analysis with Stan ----------------------------------------------------------------------------------------------



# compile stan model --------------------------------------------------------------------------------------------------
if (runStan){
  model1PL_Stan = cmdstan_model(stan_file = "Stan/01_model01.stan")

  # build stan data file ----------------------------------------------------------------------------------------------
  model1PL_Data = list(
    nObs = nrow(modelingData),
    nItems = 10,
    Y = t(correctResponseData),
    priorMeanMu = rep(0, nItems),
    priorCovMu = 10 * diag(nItems), 
    priorMeanLambda = 0,
    priorSDLambda = 10
  )

  # create starting values for theta that keep direction of theta meaningful ------------------------------------------
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = sum)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore

  # estimate model with MCMC sampling ---------------------------------------------------------------------------------
  model1PL_Samples = model1PL_Stan$sample(
    data = model1PL_Data,
    seed = 0516202301,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000,
    init = function() list(theta = startTheta)
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model1PL_Summary = model1PL_Samples$summary()
  model1PL_Draws = model1PL_Samples$draws()
  model1PL_Loo = model1PL_Samples$loo(variables = "personLike")  
  model1PL_Waic = waic(x = model1PL_Samples$draws(variables = "personLike", format = "draws_matrix"))

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples1PL = model1PL_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc1PL = buildTetrachoricPPMC(predictiveSamples1PL, rawData)

  save(model1PL_Summary, model1PL_Draws, model1PL_Loo, ppmc1PL,
       model1PL_Waic, file = "models/01_model01.RData")

} else {
  load("models/01_model01.RData")
}

# checking convergence
max(model1PL_Summary$rhat[grep(x = model1PL_Summary$variable, pattern = "mu|theta|lambda")])
 
# item parameter results
mcmc_trace(model1PL_Draws, pars = "lambda")
mcmc_dens(model1PL_Draws, pars = "lambda")

# item parameter results
print(model1PL_Summary[grep(pattern = "mu", x = model1PL_Summary$variable),], n = Inf)
print(model1PL_Summary[grep(pattern = "lambda", x = model1PL_Summary$variable),], n = Inf)

# show model fit results
View(ppmc1PL$PPMCstats)
heatmap(ppmc1PL$residTTC)

# estimate 2PL
if (runStan){
  
  model2PL_Stan = cmdstan_model(stan_file = "Stan/01_model02.stan")

  # build stan data file
  model2PL_Data = list(
    nObs = nrow(modelingData),
    nItems = 10,
    Y = t(correctResponseData),
    priorMeanMu = rep(0, nItems),
    priorCovMu = 10 * diag(nItems),
    priorMeanLambda = rep(0, nItems),
    priorCovLambda = 10 * diag(nItems)
  )

  # create starting values for theta that keep direction of theta meaningful ------------------------------------------
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = sum)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore

  # estimate 2PL model with MCMC
  model2PL_Samples = model2PL_Stan$sample(
    data = model2PL_Data,
    seed = 0517202302,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000,
    init = function() list(theta = startTheta)
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model2PL_Summary = model2PL_Samples$summary()
  model2PL_Draws = model2PL_Samples$draws()
  model2PL_Loo = model2PL_Samples$loo(variables = "personLike")  
  model2PL_Waic = waic(x = model2PL_Samples$draws(variables = "personLike", format = "draws_matrix"))

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples2PL = model2PL_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc2PL = buildTetrachoricPPMC(predictiveSamples2PL, rawData)

  save(model2PL_Summary, model2PL_Draws, model2PL_Loo, ppmc2PL, model2PL_Waic, file = "models/01_model02.RData")

} else {
  load("models/01_model02.RData")
}

# checking convergence
max(model2PL_Summary$rhat[grep(x = model2PL_Summary$variable, pattern = "mu|theta|lambda")], na.rm = TRUE)

# item parameter results
mcmc_trace(model2PL_Draws, pars = paste0("lambda[", 1:10, "]"))
mcmc_trace(model2PL_Draws, pars = paste0("mu[", 1:10, "]"))

# item parameter results
print(model2PL_Summary[grep(pattern = "mu", x = model2PL_Summary$variable),], n = Inf)
print(model2PL_Summary[grep(pattern = "lambda", x = model2PL_Summary$variable),], n = Inf)

# show model fit results
View(ppmc2PL$PPMCstats)
heatmap(ppmc2PL$residTTC)

loo_compare(list(model1PL = model1PL_Loo, 
                 model2PL = model2PL_Loo))

# model comparisons with WAIC
model1PL_Waic
model2PL_Waic

# estimate 2PL with marker item --------------------------------------------------------------------------------------
if (runStan){
  
  model2PLmarker_Stan = cmdstan_model(stan_file = "Stan/01_model03.stan")

  # build stan data file
  model2PLmarker_Data = list(
    nObs = nrow(modelingData),
    nItems = 10,
    Y = t(correctResponseData),
    priorMeanMu = rep(0, nItems),
    priorCovMu = 10 * diag(nItems),
    priorMeanLambda = rep(0, nItems-1),
    priorCovLambda = 10 * diag(nItems-1),
    priorThetaSDmean = 0,
    priorThetaSDsd = 1
  )

  # create starting values for theta that keep direction of theta meaningful ------------------------------------------
  startTheta = apply(correctResponseData, MARGIN = 1, FUN = sum)
  meanScore = mean(startTheta)
  sdScore = sd(startTheta)
  startTheta = (startTheta - meanScore) / sdScore

  # estimate 2PL model with MCMC
  model2PLmarker_Samples = model2PLmarker_Stan$sample(
    data = model2PLmarker_Data,
    seed = 0517202303,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 1000,
    init = function() list(theta = startTheta)
  )


  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  model2PLmarker_Summary = model2PLmarker_Samples$summary()
  model2PLmarker_Draws = model2PLmarker_Samples$draws()
  model2PLmarker_Loo = model2PLmarker_Samples$loo(variables = "personLike")  
  model2PLmarker_Waic = waic(x = model2PLmarker_Samples$draws(variables = "personLike", format = "draws_matrix"))

  # build model fit object --------------------------------------------------------------------------------------------
  rawData = correctResponseData
  predictiveSamples2PLmarker = model2PLmarker_Samples$draws(variables = "simY", format = "draws_matrix")

  ppmc2PLmarker = buildTetrachoricPPMC(predictiveSamples2PLmarker, rawData)

  save(
    model2PLmarker_Summary, 
    model2PLmarker_Draws, 
    model2PLmarker_Loo, 
    ppmc2PLmarker, 
    file = "models/01_model03.RData"
  )

} else {
  load("models/01_model03.RData")
}

# checking convergence
max(model2PLmarker_Summary$rhat[grep(x = model2PLmarker_Summary$variable, pattern = "mu|theta|lambda")], na.rm = TRUE)

# item parameter results
mcmc_trace(model2PLmarker_Draws, pars = paste0("lambda[", 1:9, "]"))
mcmc_trace(model2PLmarker_Draws, pars = paste0("mu[", 1:10, "]"))
mcmc_trace(model2PLmarker_Draws, pars = "thetaSD")

# item parameter results
print(model2PLmarker_Summary[grep(pattern = "mu|lambda|thetaSD", x = model2PLmarker_Summary$variable),], n = Inf)

# show model fit results
View(ppmc2PLmarker$PPMCstats)
heatmap(ppmc2PLmarker$residTTC)

loo_compare(list(model1PL = model1PL_Loo, 
                 model2PL = model2PL_Loo,
                 model2PLmarker = model2PLmarker_Loo))

# model comparisons with WAIC
model1PL_Waic
model2PL_Waic
model2PLmarker_Waic

