data {
  int<lower=0> N;           // total number of observations
  int<lower=0> nSchools;       // number of unique level-2 units in data
  array[N] int obsSchoolID;    // the level-2 ID corresponding to each row of the data
  int<lower=0> P;           // number of predictors (plus column for intercept)
  matrix[N, P] X;           // model.matrix() from R 
  vector[N] y;              // outcome
  vector[P] priorMeanBeta;       // prior mean vector for coefficients
  matrix[P, P] priorCovBeta;     // prior covariance matrix for coefficients
  real priorSigmaMean;         // prior mean parameter for residual standard deviation
  real priorSigmaSD;           // prior sd parameter for residual standard deviation
  int<lower=0> nContrasts;         // number of contrasts to estimate
  matrix[nContrasts,P] contrastMatrix;   // contrast matrix for additional effects
  int randomSlopeColumn;             // column of X that contains the random slope
  vector[2] priorTauMean;                // prior mean parameter for random effects standard deviation
  vector[2] priorTauSD;                  // prior sd parameter for random effects standard deviation
  real priorRandomEffectsCorrLJK;       // prior for random effects correlation  
}

parameters {
  vector[P] beta;                    // vector of coefficients for beta
  real<lower=0> sigma;               // residual standard deviation
  array[nSchools] vector[2] randomEffects;     // random intercept and slope for each level-2 unit
  cholesky_factor_corr[2] randomEffectsCorrL;
  vector<lower=0>[2] randomEffectsSD;
}

model {
  // for random effects;
  vector[2] meanRandomEffects = rep_vector(0, 2);
  matrix[2, 2] randomEffectsCovL;
  randomEffectsCorrL ~ lkj_corr_cholesky(priorRandomEffectsCorrLJK);
  randomEffectsSD ~ lognormal(priorTauMean, priorTauSD);
  randomEffectsCovL = diag_pre_multiply(randomEffectsSD, randomEffectsCorrL);

  beta ~ multi_normal(priorMeanBeta, priorCovBeta);        // prior for coefficients
  sigma ~ lognormal(priorSigmaMean, priorSigmaSD);         // prior for sigma

  for (school in 1:nSchools){
    randomEffects[school] ~ multi_normal_cholesky(meanRandomEffects, randomEffectsCovL); 
  }

  for (obs in 1:N){
    y[obs] ~ normal(
      X[obs,]*beta + randomEffects[obsSchoolID[obs],1] + randomEffects[obsSchoolID[obs],2]*X[obs,randomSlopeColumn], 
      sigma); 
  }
  
}

generated quantities {

  vector[N] personLike = rep_vector(0.0, N);
  for (obs in 1:N){
    personLike[obs] = 
      normal_lpdf(y[obs] | 
                  X[obs,]*beta + randomEffects[obsSchoolID[obs],1] + randomEffects[obsSchoolID[obs],2]*X[obs,randomSlopeColumn],
                  sigma);
  }

  // contrast estimation (linear combinations)
  vector[nContrasts] constrastEstimates;
  constrastEstimates = contrastMatrix*beta;

  // transform correlation matrix and SD vector to covariance matrix
  corr_matrix[2] randomEffectsCorr;
  cov_matrix[2] randomEffectsCov; 
  matrix[2, 2] randomEffectsCovL;
 
  randomEffectsCorr = multiply_lower_tri_self_transpose(randomEffectsCorrL);
  randomEffectsCovL = diag_pre_multiply(randomEffectsSD, randomEffectsCorrL);
  randomEffectsCov = multiply_lower_tri_self_transpose(randomEffectsCovL);
}
