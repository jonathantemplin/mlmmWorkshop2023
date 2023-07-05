data {
  int<lower=0> N;         // total number of observations
  int<lower=0> P;         // number of predictors (plus column for intercept)
  matrix[N, P] X;         // model.matrix() from R 
  vector[N] y;            // outcome
  vector[P] priorMeanBeta;     // prior mean vector for coefficients
  matrix[P, P] priorCovBeta;   // prior covariance matrix for coefficients
  real priorSigmaMean;         // prior rate parameter for residual standard deviation
  real priorSigmaSD;           // prior shape parameter for residual standard deviation
}

parameters {
  vector[P] beta;         // vector of coefficients for beta
  real<lower=0> sigma;    // residual standard deviation
}

model {
  beta ~ multi_normal(priorMeanBeta, priorCovBeta); // prior for coefficients
  sigma ~ lognormal(priorSigmaMean, priorSigmaSD);         // prior for sigma
  y ~ normal(X*beta, sigma);              // linear model predicting single outcome
}

generated quantities {

  vector[N] personLike = rep_vector(0.0, N);

  for (i in 1:N) {
    personLike[i] = normal_lpdf(y[i] | X[i]*beta, sigma);
  }

}
