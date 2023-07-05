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
  real priorTauMean;                // prior mean parameter for random intercept standard deviation
  real priorTauSD;                  // prior sd parameter for random intercept standard deviation
  int<lower=0> nContrasts;         // number of contrasts to estimate
  matrix[nContrasts,P] contrastMatrix;   // contrast matrix for additional effects
}

parameters {
  vector[P] beta;                    // vector of coefficients for beta
  real<lower=0> sigma;               // residual standard deviation
  vector[nSchools] randomIntercept;     // random intercept for each level-2 unit
  real<lower=0> tauIntercept;           // random intercept for each level-2 unit
}

model {
  beta ~ multi_normal(priorMeanBeta, priorCovBeta);        // prior for coefficients
  sigma ~ lognormal(priorSigmaMean, priorSigmaSD);         // prior for sigma
  tauIntercept ~ lognormal(priorTauMean, priorTauSD);      // prior for tau

  for (school in 1:nSchools){
    randomIntercept[school] ~ normal(0, tauIntercept);
  }

  for (obs in 1:N){
    y[obs] ~ normal(X[obs,]*beta + randomIntercept[obsSchoolID[obs]], sigma); // linear model
  }
  
}

generated quantities {
  real ICC;
  ICC = tauIntercept^2/(tauIntercept^2 + sigma^2);

  vector[N] personLike = rep_vector(0.0, N);
  for (obs in 1:N){
    personLike[obs] = 
      normal_lpdf(y[obs] | X[obs,]*beta + randomIntercept[obsSchoolID[obs]], sigma);
  }

  // contrast estimation (linear combinations)
  vector[nContrasts] constrastEstimates;
  constrastEstimates = contrastMatrix*beta;

}
