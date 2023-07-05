data {
  int<lower=0> N;         // total number of observations
  int<lower=0> P;         // number of predictors (plus column for intercept)
  matrix[N, P] X;         // model.matrix() from R 
  array[N] int y;         // outcome (now an array instead of vector)
  vector[P] priorMeanBeta;     // prior mean vector for coefficients
  matrix[P, P] priorCovBeta;   // prior covariance matrix for coefficients
}

parameters {
  vector[P] beta;         // vector of coefficients for beta
}

model {
  beta ~ multi_normal(priorMeanBeta, priorCovBeta); // prior for coefficients
  y ~ bernoulli_logit(X*beta); // linear model predicting single outcome
}

generated quantities {
  real prob;
  prob = exp(beta[1])/(1+exp(beta[1]));

  vector[N] personLike = rep_vector(0.0, N);  
  for (i in 1:N) {
    personLike[i] = bernoulli_logit_lpmf(y[i] | X[i]*beta);
  }
}
