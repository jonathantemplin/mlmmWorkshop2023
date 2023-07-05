# Multilevel Modeling: Predicting Observed Sum Score for Students Nested in Schools 
# Using Bayesian Estimation 

# clear workspace =====================================================================================================
rm(list = ls())

# set options =========================================================================================================
# Set width of output and number of significant digits printed,
# number of digits before using scientific notation, shut off significance stars, # cores
options(width = 120, digits = 8, scipen = 9, show.signif.stars = FALSE, mc.cores = 4)

runStan = TRUE # set to TRUE to run Stan models; if FALSE, will load saved results


# install functions needed for analysis ===============================================================================
addUnitMeans = function(data, unitVariable, meanVariables, newNames){
  # create unit-level data set:

  unitMeans = t(sapply(
    X = unique(data[,unitVariable]),
    FUN = function(x, data) {
      return(c(
        x,
        length(which(data[,unitVariable] == x)),
        apply(
          X = as.data.frame(data[which(data[,unitVariable] == x), meanVariables]), 
          MARGIN = 2, 
          FUN = mean, 
          rm.na=TRUE
          )
        )
      )
    },data = data
  ))

  unitMeans = data.frame(unitMeans)
  names(unitMeans) = c(unitVariable, paste0("Nper", unitVariable), newNames)

  newData = merge(x = data, y = unitMeans, by = unitVariable)
  return(newData)
}


# Package installation ================================================================================================

needed_packages = c("psych", "ggplot2", "cmdstanr", "HDInterval", "bayesplot", "loo")
for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage == FALSE) {
    install.packages(needed_packages[i])
  }
  library(needed_packages[i], character.only = TRUE)
}


# =====================================================================================================================
# Data Import, Manipulation, and Description

# Load example data from folder within working directory
load("data/modelingData.RData")

# Create school means for sum score student outcome and free/reduced lunch student predictor 
modelingData = addUnitMeans(data = modelingData, unitVariable = "schoolID", 
                            meanVariables = c("sumScore", "frlunch"), 
                            newNames = c("SMsumScore", "SMfrlunch"))


# Descriptive statistics for student variables and new school means
print(
  describe(
    modelingData[c("sumScore", "frlunch", "SMsumScore", "SMfrlunch")]
    ),
  digits = 3
)

# Center school lunch near mean to use as observed level-2 predictor
modelingData$SMfrlunch30 = modelingData$SMfrlunch - .30

# Cluster-mean-center student lunch at school mean to use as observed level-1 predictor
modelingData$WSfrlunch = modelingData$frlunch - modelingData$SMfrlunch

# Preliminary data steps for two-level models in Stan
# Create new unique school IDs starting at 1
schoolIDs = unique(modelingData$schoolID)
nSchools = length(schoolIDs)
# Build Stan's schoolID ranging from 1:nSchools
stanSchoolIDs = 1:nSchools
schoolIDtable = data.frame(schoolID = schoolIDs, stanSchoolID = stanSchoolIDs)
# Merge into data
modelingData = merge(x = modelingData, y = schoolIDtable, by = "schoolID")


# =====================================================================================================================
# Partitioning Variance in the Sum Score Outcome using General Models
  # Single-level empty model predicting observed sum score ignoring school
  # stan syntax in 02_model01.stan


if (runStan){
  modelEmptyGLM_Stan = cmdstan_model(stan_file = "Stan/02_model01.stan")

  # Write model for the means formula
  modelEmptyGLM_Formula = formula(sumScore ~ 1, data = modelingData)

  # Make model matrix
  modelEmptyGLM_modelMatrix = model.matrix(modelEmptyGLM_Formula, data = modelingData)
  dim(modelEmptyGLM_modelMatrix)

  # Build Stan data from model matrix
  modelEmptyGLM_Data = list(
      X = modelEmptyGLM_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelEmptyGLM_modelMatrix), 
      P = ncol(modelEmptyGLM_modelMatrix), 
      priorMeanBeta = rep(0, ncol(modelEmptyGLM_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelEmptyGLM_modelMatrix), 
                          ncol = ncol(modelEmptyGLM_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100
  )

  # Run Stan file
  modelEmptyGLM_Samples = modelEmptyGLM_Stan$sample(
      data = modelEmptyGLM_Data, 
      seed = 0608202301, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelEmptyGLM_Summary = modelEmptyGLM_Samples$summary()
  modelEmptyGLM_Draws = modelEmptyGLM_Samples$draws()
  modelEmptyGLM_Loo = modelEmptyGLM_Samples$loo(variables = "personLike")  

  save(modelEmptyGLM_Summary, modelEmptyGLM_Draws, modelEmptyGLM_Loo, file = "models/02_model01.RData")

} else {
  load("models/02_model01.RData")
}

# Assess convergence: summary of all parameters
max(modelEmptyGLM_Summary$rhat[grep(pattern = "beta|sigma", x = modelEmptyGLM_Summary$variable)])

# parameter results
print(modelEmptyGLM_Summary[grep(pattern = "beta|sigma", x = modelEmptyGLM_Summary$variable),], n = Inf)


# =====================================================================================================================
# Two-level empty model predicting observed sum score with students nested in schools

if (runStan){

  # Making it easier to supply input into Stan (new parts with indented comments below)
  modelRandomIntercept_Stan = cmdstan_model(stan_file = "Stan/02_model02.stan")

  # Write model for the means formula
  modelRandomIntercept_Formula = formula(sumScore ~ 1, data = modelingData)

  # Make model matrix
  modelRandomIntercept_modelMatrix = model.matrix(modelRandomIntercept_Formula, data = modelingData)


  # Build Stan data from model matrix
  modelRandomIntercept_Data = list(
      X = modelRandomIntercept_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelRandomIntercept_modelMatrix),
      P = ncol(modelRandomIntercept_modelMatrix), 
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      priorMeanBeta = rep(0, ncol(modelRandomIntercept_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelRandomIntercept_modelMatrix), 
                          ncol = ncol(modelRandomIntercept_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      priorTauMean = 0,
      priorTauSD = 100,
      nContrasts = 0,
      contrastMatrix = matrix(data = NA, nrow = 0, ncol = ncol(modelRandomIntercept_modelMatrix))
  )

  # Run Stan file
  modelRandomIntercept_Samples = modelRandomIntercept_Stan$sample(
      data = modelRandomIntercept_Data, 
      seed = 0608202302, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelRandomIntercept_Summary = modelRandomIntercept_Samples$summary()
  modelRandomIntercept_Draws = modelRandomIntercept_Samples$draws()
  modelRandomIntercept_Loo = modelRandomIntercept_Samples$loo(variables = "personLike")

  save(modelRandomIntercept_Summary, modelRandomIntercept_Draws, modelRandomIntercept_Loo, 
       file = "models/02_model02a.RData")
} else {
  load("models/02_model02a.RData")
}

# Assess convergence: summary of all parameters
max(
  modelRandomIntercept_Summary$rhat[
    grep(pattern = "beta|sigma|tauIntercept", x = modelRandomIntercept_Summary$variable)
    ]
)

# parameter results
print(
  modelRandomIntercept_Summary[
    grep(pattern = "beta|sigma|tauIntercept|ICC|contrastEstimates", x = modelRandomIntercept_Summary$variable),
    ], 
  n = Inf
)

mcmc_dens(modelRandomIntercept_Draws, pars = "tauIntercept")

# ICC Results
mcmc_trace(modelRandomIntercept_Draws, pars = "ICC")
mcmc_dens(modelRandomIntercept_Draws, pars = "ICC")

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLM_Loo, 
    emptyRandomIntercept = modelRandomIntercept_Loo
  )
)

# =====================================================================================================================
# Partitioning Variance in the Binary Free/Reduced Lunch Predictor using Generalized Models

# Single-level empty model predicting observed free/reduced lunch ignoring school

if (runStan){
  modelEmptyGLMfr_Stan = cmdstan_model(stan_file = "Stan/02_model03.stan")

  # Write model for the means formula
  modelEmptyGLMfr_Formula = formula(frlunch ~ 1, data = modelingData)

  # Make model matrix
  modelEmptyGLMfr_modelMatrix = model.matrix(modelEmptyGLMfr_Formula, data = modelingData)


  # Build Stan data from model matrix
  modelEmptyGLMfr_Data = list(
    X = modelEmptyGLMfr_modelMatrix, 
    y = modelingData$frlunch,
    N = nrow(modelEmptyGLMfr_modelMatrix), 
    P = ncol(modelEmptyGLMfr_modelMatrix), 
    priorMeanBeta = rep(0, ncol(modelEmptyGLMfr_modelMatrix)), 
    priorCovBeta = diag(
      x = 10000, 
      nrow = ncol(modelEmptyGLMfr_modelMatrix), 
      ncol = ncol(modelEmptyGLMfr_modelMatrix)
      )
  )

  # Run Stan file
  modelEmptyGLMfr_Samples = modelEmptyGLMfr_Stan$sample(
    data = modelEmptyGLMfr_Data, 
    seed = 0608202304, 
    chains = 4, 
    parallel_chains = 4, 
    iter_warmup = 2000, 
    iter_sampling = 1000,
    refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelEmptyGLMfr_Summary = modelEmptyGLMfr_Samples$summary()
  modelEmptyGLMfr_Draws = modelEmptyGLMfr_Samples$draws()
  modelEmptyGLMfr_Loo = modelEmptyGLMfr_Samples$loo(variables = "personLike")

  save(modelEmptyGLMfr_Summary, modelEmptyGLMfr_Draws, modelEmptyGLMfr_Loo, 
       file = "models/02_model03.RData")

} else {
  load("models/02_model03.RData")
}


# Assess convergence: summary of all parameters
max(
  modelEmptyGLMfr_Summary$rhat[
    grep(pattern = "beta", x = modelEmptyGLMfr_Summary$variable)
    ]
)

# parameter results
print(
  modelEmptyGLMfr_Summary[
    grep(pattern = "beta|prob", x = modelEmptyGLMfr_Summary$variable),
    ], 
  n = Inf
)


# =====================================================================================================================
# Two-level empty model predicting observed free/reduced lunch with students nested in schools

if (runStan){
    
  modelEmptyGLMfrRandomIntercept_Stan = cmdstan_model(stan_file = "Stan/02_model04.stan")

  # Write model for the means formula
  modelEmptyGLMfrRandomIntercept_Formula = formula(frlunch ~ 1, data = modelingData)

  # Make model matrix
  modelEmptyGLMfrRandomIntercept_modelMatrix = 
    model.matrix(modelEmptyGLMfrRandomIntercept_Formula, data = modelingData)

  # Build Stan data from model matrix
  modelEmptyGLMfrRandomIntercept_Data = list(
    X = modelEmptyGLMfrRandomIntercept_modelMatrix, 
    y = modelingData$frlunch,
    N = nrow(modelingData), 
    P = ncol(modelEmptyGLMfrRandomIntercept_modelMatrix), 
    nSchools = nSchools, 
    obsSchoolID = modelingData$stanSchoolID,
    priorMeanBeta = rep(0, ncol(modelEmptyGLMfrRandomIntercept_modelMatrix)), 
    priorCovBeta = diag(
      x = 10000, 
      nrow = ncol(modelEmptyGLMfrRandomIntercept_modelMatrix), 
      ncol = ncol(modelEmptyGLMfrRandomIntercept_modelMatrix)
      ),
    priorTauMean = 0,
    priorTauSD = 100
  )

  # Run Stan file
  modelEmptyGLMfrRandomIntercept_Samples = 
    modelEmptyGLMfrRandomIntercept_Stan$sample(
      data = modelEmptyGLMfrRandomIntercept_Data, 
      seed = 0608202305, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
    )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelEmptyGLMfrRandomIntercept_Summary = modelEmptyGLMfrRandomIntercept_Samples$summary()
  modelEmptyGLMfrRandomIntercept_Draws = modelEmptyGLMfrRandomIntercept_Samples$draws()
  modelEmptyGLMfrRandomIntercept_Loo = modelEmptyGLMfrRandomIntercept_Samples$loo(variables = "personLike")

  save(modelEmptyGLMfrRandomIntercept_Summary, modelEmptyGLMfrRandomIntercept_Draws, 
       modelEmptyGLMfrRandomIntercept_Loo, file = "models/02_model04.RData")

} else {
  load("models/02_model04.RData")
}

# Assess convergence: summary of all parameters
max(
  modelEmptyGLMfrRandomIntercept_Summary$rhat[
    grep(pattern = "beta|tauIntercept", x = modelEmptyGLMfrRandomIntercept_Summary$variable)
    ]
)

# parameter results
print(
  modelEmptyGLMfrRandomIntercept_Summary[
    grep(pattern = "beta|tauIntercept|ICC|prob", x = modelEmptyGLMfrRandomIntercept_Summary$variable),
    ], 
  n = Inf
)

# ICC Results
mcmc_trace(modelEmptyGLMfrRandomIntercept_Draws, pars = "ICC")
mcmc_dens(modelEmptyGLMfrRandomIntercept_Draws, pars = "ICC")

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLMfr_Loo, 
    emptyRandomIntercept = modelEmptyGLMfrRandomIntercept_Loo
  )
)

#########################################################################################################
# Smushed model predicting observed sum score with students nested in schools
# uses same Stan model as 02_model02.stan: Random intercepts for schools

if (runStan){

  # Write model for themeans formula
  modelSmushed_Formula = formula(sumScore ~ 1 + frlunch, data = modelingData)

  # Make model matrix
  modelSmushed_modelMatrix = model.matrix(modelSmushed_Formula, data = modelingData)
colnames(modelSmushed_modelMatrix)
  # Build Stan data from model matrix
  modelSmushed_Data = list(
      X = modelSmushed_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelSmushed_modelMatrix),
      P = ncol(modelSmushed_modelMatrix), 
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      priorMeanBeta = rep(0, ncol(modelSmushed_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelSmushed_modelMatrix), 
                          ncol = ncol(modelSmushed_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      priorTauMean = 0,
      priorTauSD = 100,
      nContrasts = 0,
      contrastMatrix = matrix(data = NA, nrow = 0, ncol = ncol(modelSmushed_modelMatrix))
  )

  # Run Stan file
  modelSmushed_Samples = modelRandomIntercept_Stan$sample(
      data = modelSmushed_Data, 
      seed = 0608202303, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelSmushed_Summary = modelSmushed_Samples$summary()
  modelSmushed_Draws = modelSmushed_Samples$draws()
  modelSmushed_Loo = modelSmushed_Samples$loo(variables = "personLike")

  save(modelSmushed_Summary, modelSmushed_Draws, modelSmushed_Loo, file = "models/02_model02b.RData")

} else {
  load("models/02_model02b.RData")
}

# Assess convergence: summary of all parameters
max(
  modelSmushed_Summary$rhat[
    grep(pattern = "beta|sigma|tauIntercept", x = modelSmushed_Summary$variable)
    ]
)

# parameter results
print(
  modelSmushed_Summary[
    grep(pattern = "beta|sigma|tauIntercept|ICC", x = modelSmushed_Summary$variable),
    ], 
  n = Inf
)

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLM_Loo, 
    emptyRandomIntercept = modelRandomIntercept_Loo,
    smushedRIfrLunch = modelSmushed_Loo
  )
)

# =========================================================================================================
# Add centered school mean to unsmush the level-1 slope for frlunch
# uses same Stan model as 02_model03.stan: Random intercepts for schools

if (runStan){

  # Write model for the means formula
  modelLunch_Formula = formula(sumScore ~ 1 + frlunch + SMfrlunch30, data = modelingData)

  # Make model matrix
  modelLunch_modelMatrix = model.matrix(modelLunch_Formula, data = modelingData)

  modelLunch_contrastMatrix = matrix(
      data = c(0, 1, 1), 
      nrow = 1
  )

  # Build Stan data from model matrix
  modelLunch_Data = list(
      X = modelLunch_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelLunch_modelMatrix),
      P = ncol(modelLunch_modelMatrix), 
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      priorMeanBeta = rep(0, ncol(modelLunch_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelLunch_modelMatrix), 
                          ncol = ncol(modelLunch_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      priorTauMean = 0,
      priorTauSD = 100,
      nContrasts = nrow(modelLunch_contrastMatrix),
      contrastMatrix = modelLunch_contrastMatrix
  )

  # Run Stan file
  modelLunch_Samples = modelRandomIntercept_Stan$sample(
      data = modelLunch_Data, 
      seed = 0608202306, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelLunch_Summary = modelLunch_Samples$summary()
  modelLunch_Draws = modelLunch_Samples$draws()
  modelLunch_Loo = modelLunch_Samples$loo(variables = "personLike")

  save(modelLunch_Summary, modelLunch_Draws, modelLunch_Loo, file = "models/02_model02c.RData")

} else {
  load("models/02_model02c.RData")
}

# Assess convergence: summary of all parameters
max(
  modelLunch_Summary$rhat[
    grep(pattern = "beta|sigma|tauIntercept", x = modelLunch_Summary$variable)
    ]
)

# parameter results
print(
  modelLunch_Summary[
    grep(pattern = "beta|sigma|tauIntercept|ICC|constrastEstimates", x = modelLunch_Summary$variable),
    ], 
  n = Inf
)

mcmc_dens(modelLunch_Draws, pars = "constrastEstimates[1]")

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLM_Loo, 
    emptyRandomIntercept = modelRandomIntercept_Loo,
    smushedRIfrLunch = modelSmushed_Loo,
    correctRIfrLunch = modelLunch_Loo
  )
)

# =========================================================================================================
# Cluster-mean-centered version of modelLunch
# uses same Stan model as 02_model03.stan: Random intercepts for schools

if (runStan){

  # Write model for the means formula
  modelLunchCMC_Formula = formula(sumScore ~ 1 + WSfrlunch + SMfrlunch30, data = modelingData)

  # Make model matrix
  modelLunchCMC_modelMatrix = model.matrix(modelLunchCMC_Formula, data = modelingData)

  modelLunchCMC_contrastMatrix = matrix(
      data = c(0, -1, 1), 
      nrow = 1
  )

  # Build Stan data from model matrix
  modelLunchCMC_Data = list(
      X = modelLunchCMC_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelLunchCMC_modelMatrix),
      P = ncol(modelLunchCMC_modelMatrix), 
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      priorMeanBeta = rep(0, ncol(modelLunchCMC_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelLunchCMC_modelMatrix), 
                          ncol = ncol(modelLunchCMC_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      priorTauMean = 0,
      priorTauSD = 100,
      nContrasts = nrow(modelLunchCMC_contrastMatrix),
      contrastMatrix = modelLunchCMC_contrastMatrix
  )

  # Run Stan file
  modelLunchCMC_Samples = modelRandomIntercept_Stan$sample(
      data = modelLunchCMC_Data, 
      seed = 0608202307, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelLunchCMC_Summary = modelLunchCMC_Samples$summary()
  modelLunchCMC_Draws = modelLunchCMC_Samples$draws()
  modelLunchCMC_Loo = modelLunchCMC_Samples$loo(variables = "personLike")

  save(modelLunchCMC_Summary, modelLunchCMC_Draws, modelLunchCMC_Loo, file = "models/02_model02d.RData")

} else {
  load("models/02_model02d.RData")
}

# Assess convergence: summary of all parameters
max(
  modelLunchCMC_Summary$rhat[
    grep(pattern = "beta|sigma|tauIntercept", x = modelLunchCMC_Summary$variable)
    ]
)

# parameter results
print(
  modelLunchCMC_Summary[
    grep(pattern = "beta|sigma|tauIntercept|ICC|constrastEstimates", x = modelLunchCMC_Summary$variable),
    ], 
  n = Inf
)

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLM_Loo, 
    emptyRandomIntercept = modelRandomIntercept_Loo,
    smushedRIfrLunch = modelSmushed_Loo,
    correctRIfrLunch = modelLunch_Loo,
    CMCRIfrLunch = modelLunchCMC_Loo
  )
)

# =========================================================================================================
#  Models with a Random Slope across Schools for Cluster-Mean-Centered Student Free/Reduced Lunch 

if (runStan){

  modelRandomSlope_Stan = cmdstan_model(stan_file = "Stan/02_model05.stan")

# Write model for the means formula
  modelRandomSlope_Formula = formula(sumScore ~ 1 + WSfrlunch + SMfrlunch30, data = modelingData)

  # Make model matrix
  modelRandomSlope_modelMatrix = model.matrix(modelRandomSlope_Formula, data = modelingData)

  # Build Stan data from model matrix
  modelRandomSlope_Data = list(
      X = modelRandomSlope_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelRandomSlope_modelMatrix),
      P = ncol(modelRandomSlope_modelMatrix), 
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      priorMeanBeta = rep(0, ncol(modelRandomSlope_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelRandomSlope_modelMatrix), 
                          ncol = ncol(modelRandomSlope_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      nContrasts = 0,
      contrastMatrix = matrix(data = NA, nrow = 0, ncol = ncol(modelRandomSlope_modelMatrix)),
      randomSlopeColumn = which(colnames(modelRandomSlope_modelMatrix) == "WSfrlunch"),
      priorTauMean = c(0, -5),
      priorTauSD = c(100, .5),
      priorRandomEffectsCorrLJK = 1
  )

  # Run Stan file
  modelRandomSlope_Samples = modelRandomSlope_Stan$sample(
      data = modelRandomSlope_Data, 
      seed = 0608202308, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelRandomSlope_Summary = modelRandomSlope_Samples$summary()
  modelRandomSlope_Draws = modelRandomSlope_Samples$draws()
  modelRandomSlope_Loo = modelRandomSlope_Samples$loo(variables = "personLike")

  save(modelRandomSlope_Summary, modelRandomSlope_Draws, modelRandomSlope_Loo, file = "models/02_model05a.RData")

} else {
  load("models/02_model05a.RData")
}

# Assess convergence: summary of all parameters
max(
  modelRandomSlope_Summary$rhat[
    grep(pattern = "beta|sigma|randomEffectsSD|randomEffectsCorr", x = modelRandomSlope_Summary$variable)
    ],
  na.rm = TRUE
)


# parameter results
print(
  modelRandomSlope_Summary[
    grep(pattern = "beta|sigma|randomEffectsSD|randomEffectsCorr|randomEffectsCov", x = modelRandomSlope_Summary$variable),
    ],
  n = Inf
)

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLM_Loo, 
    emptyRandomIntercept = modelRandomIntercept_Loo,
    smushedRIfrLunch = modelSmushed_Loo,
    correctRIfrLunch = modelLunch_Loo,
    CMCRIfrLunch = modelLunchCMC_Loo,
    randomSlope = modelRandomSlope_Loo
  )
)

# =====================================================================================================================
#  Models with a Random Slope across Schools for Cluster-Mean-Centered Student Free/Reduced Lunch with 
#  Cross-level interaction
#  uses same Stan model as 02_model05.stan: Random intercepts and slope for schools


if (runStan){

# Write model for the means formula
  modelRandomSlopeCLI_Formula = 
    formula(sumScore ~ 1 + WSfrlunch + SMfrlunch30 + WSfrlunch:SMfrlunch30, data = modelingData)

  # Make model matrix
  modelRandomSlopeCLI_modelMatrix = model.matrix(modelRandomSlopeCLI_Formula, data = modelingData)

  # Build Stan data from model matrix
  modelRandomSlopeCLI_Data = list(
      X = modelRandomSlopeCLI_modelMatrix, 
      y = modelingData$sumScore,
      N = nrow(modelRandomSlopeCLI_modelMatrix),
      P = ncol(modelRandomSlopeCLI_modelMatrix), 
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      priorMeanBeta = rep(0, ncol(modelRandomSlopeCLI_modelMatrix)), 
      priorCovBeta = diag(x = 10000, 
                          nrow = ncol(modelRandomSlopeCLI_modelMatrix), 
                          ncol = ncol(modelRandomSlopeCLI_modelMatrix)), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      nContrasts = 0,
      contrastMatrix = matrix(data = NA, nrow = 0, ncol = ncol(modelRandomSlopeCLI_modelMatrix)),
      randomSlopeColumn = which(colnames(modelRandomSlopeCLI_modelMatrix) == "WSfrlunch"),
      priorTauMean = c(0, -5),
      priorTauSD = c(100, .5),
      priorRandomEffectsCorrLJK = 1
  )

  # Run Stan file
  modelRandomSlopeCLI_Samples = modelRandomSlope_Stan$sample(
      data = modelRandomSlopeCLI_Data, 
      seed = 0608202309, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelRandomSlopeCLI_Summary = modelRandomSlopeCLI_Samples$summary()
  modelRandomSlopeCLI_Draws = modelRandomSlopeCLI_Samples$draws()
  modelRandomSlopeCLI_Loo = modelRandomSlopeCLI_Samples$loo(variables = "personLike")

  save(modelRandomSlopeCLI_Summary, modelRandomSlopeCLI_Draws, modelRandomSlopeCLI_Loo, 
       file = "models/02_model05b.RData")

} else {
  load("models/02_model05b.RData")
}

max(
  modelRandomSlopeCLI_Summary$rhat[
    grep(pattern = "beta|sigma|randomEffectsSD|randomEffectsCorr", x = modelRandomSlopeCLI_Summary$variable)
    ],
  na.rm = TRUE
)

# parameter results
print(
  modelRandomSlopeCLI_Summary[
    grep(pattern = "beta|sigma|randomEffectsSD|randomEffectsCorr|randomEffectsCov", x = modelRandomSlopeCLI_Summary$variable),
    ],
  n = Inf
)

 # compare relative model fit
loo_compare(
  list(
    emptyGLM = modelEmptyGLM_Loo, 
    emptyRandomIntercept = modelRandomIntercept_Loo,
    smushedRIfrLunch = modelSmushed_Loo,
    correctRIfrLunch = modelLunch_Loo,
    CMCRIfrLunch = modelLunchCMC_Loo,
    randomSlope = modelRandomSlope_Loo,
    randomSlopeCLI = modelRandomSlopeCLI_Loo
  )
)
# =====================================================================================================================
# Add centered school mean to unsmush the level-1 slope for frlunch
# But--- school mean is estimated rather than set at MLE: Multivariate MLM

if (runStan){

  modelMultivariate_Stan = cmdstan_model(stan_file = "Stan/02_model06.stan")

  modelMultivariate_Data = list(
      N = nrow(modelingData),
      nSchools = nSchools, 
      obsSchoolID = modelingData$stanSchoolID,
      frlunch = modelingData$frlunch,
      score = modelingData$sumScore,
      priorBetaMean = 0, 
      priorBetaSD = sqrt(10000), 
      priorSigmaMean = 0,
      priorSigmaSD = 100,
      priorTauScoreMean = 0,
      priorTauScoreSD = 100,
      priorTauFrlunchMean = 0,
      priorTauFrlunchSD = 100
  )

  # Run Stan file
  modelMultivariate_Samples = modelMultivariate_Stan$sample(
      data = modelMultivariate_Data, 
      seed = 0608202310, 
      chains = 4, 
      parallel_chains = 4, 
      iter_warmup = 2000, 
      iter_sampling = 1000,
      refresh = 100
  )

  # save MCMC object to disk to remove time burden --------------------------------------------------------------------
  modelMultivariate_Summary = modelMultivariate_Samples$summary()
  modelMultivariate_Draws = modelMultivariate_Samples$draws()
  modelMultivariate_Loo = modelMultivariate_Samples$loo(variables = "personLike")

  save(modelMultivariate_Summary, modelMultivariate_Draws, modelMultivariate_Loo, file = "models/02_model06.RData")
  
} else {
  load("models/02_model06.RData")
}

max(
  modelMultivariate_Summary$rhat[
    grep(pattern = 
      "scoreIntercept|frlunchIntercept|frlunchSlope|frlunchMeanSlope|sigma|tauScoreIntercept|tauFrlunchIntercept|ICC",
     x = modelMultivariate_Summary$variable)
    ],
  na.rm = TRUE
)

# parameter results
print(
  modelMultivariate_Summary[
    grep(pattern = 
      "scoreIntercept|frlunchIntercept|frlunchSlope|frlunchMeanSlope|sigma|tauScoreIntercept|tauFrlunchIntercept|ICC", 
      x = modelMultivariate_Summary$variable),
    ],
  n = Inf
)
