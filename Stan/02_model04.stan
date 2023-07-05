data {
  int<lower=0> N;           // total number of observations
  int<lower=0> nSchools;       // number of unique level-2 units in data
  array[N] int obsSchoolID;    // the level-2 ID corresponding to each row of the data
  int<lower=0> P;           // number of predictors (plus column for intercept)
  matrix[N, P] X;           // model.matrix() from R 
  array[N] int y;           // outcome (now in array instead of vector)
  vector[P] priorMeanBeta;       // prior mean vector for coefficients
  matrix[P, P] priorCovBeta;     // prior covariance matrix for coefficients
  real priorTauMean;                // prior mean parameter for random intercept standard deviation
  real priorTauSD;                  // prior sd parameter for random intercept standard deviation
}

parameters {
  vector[P] beta;                    // vector of coefficients for beta
  vector[nSchools] randomIntercept;     // random intercept for each level-2 unit
  real<lower=0> tauIntercept;           // random intercept for each level-2 unit
}

model {
  beta ~ multi_normal(priorMeanBeta, priorCovBeta);         // prior for coefficients
  tauIntercept ~ lognormal(priorTauMean, priorTauSD);       // prior for random intercept standard deviation

  for (school in 1:nSchools){
    randomIntercept[school] ~ normal(0, tauIntercept);
  }

  for (obs in 1:N){
    y[obs] ~ bernoulli_logit(X[obs,]*beta + randomIntercept[obsSchoolID[obs]]); // linear model
  }

}

generated quantities {
  real prob;
  prob = exp(beta[1])/(1+exp(beta[1]));

  real ICC;
  ICC = tauIntercept^2/(tauIntercept^2 + ((pi()^2)/3));

  vector[N] personLike = rep_vector(0.0, N);

  for (n in 1:N){
    personLike[n] = bernoulli_logit_lpmf(y[n] | X[n,]*beta + randomIntercept[obsSchoolID[n]]);
  }
}
