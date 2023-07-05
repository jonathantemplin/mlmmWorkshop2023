data {
  int<lower=0> nObs;                            // number of observations
  int<lower=0> nItems;                          // number of items
  array[nItems, nObs] int<lower=0, upper=1>  Y; // item responses in an array

  vector[nItems] priorMeanMu;             // prior mean vector for intercept parameters
  matrix[nItems, nItems] priorCovMu;      // prior covariance matrix for intercept parameters

  vector[nItems-1] priorMeanLambda;         // prior mean vector for loading parameters
  matrix[nItems-1, nItems-1] priorCovLambda;   // prior covariance matrix for loading parameters

  real priorThetaSDmean;                  // prior mean for latent variable standard deviation
  real priorThetaSDsd;                    // prior sd for latent variable standard deviation
}

parameters {
  vector[nObs] theta;                // the latent variables (one for each person)
  vector[nItems] mu;                 // the item intercepts (one for each item)
  vector[nItems-1] lambda;             // the factor loadings/item discriminations (one for each item)
  real<lower=0> thetaSD;             // the standard deviation of the latent variable
}

model {

  lambda ~ multi_normal(priorMeanLambda, priorCovLambda); // Prior for item discrimination/factor loadings
  mu ~ multi_normal(priorMeanMu, priorCovMu);             // Prior for item intercepts
  thetaSD ~ lognormal(priorThetaSDmean, priorThetaSDsd);                              // Prior for latent variable standard deviation

  theta ~ normal(0, thetaSD);                         // Prior for latent variable (with mean/sd specified)
  
  Y[1] ~ bernoulli_logit(mu[1] + theta); // Likelihood for first item
  for (item in 2:nItems){
    Y[item] ~ bernoulli_logit(mu[item] + lambda[item-1]*theta);
  }
  
}


generated quantities{

  // converting item difficulties to item intercepts (easiness)
  vector[nItems] b;
  
  b[1] = -1*mu[1];
  for (item in 2:nItems){
    b[item] = -1*mu[item]/lambda[item-1];
  }
  
  // for PPMC:
  array[nItems, nObs] int<lower=0> simY;
  
  // for LOO/WAIC:
  vector[nObs] personLike = rep_vector(0.0, nObs);
  
  for (obs in 1:nObs){
    simY[1, obs] = bernoulli_logit_rng(mu[1] + theta[obs]);
    
    personLike[obs] = 
        personLike[obs] + 
        bernoulli_logit_lpmf(Y[1, obs] | mu[1] + theta[obs]);

    for (item in 2:nItems){
      // generate data based on distribution and model
      simY[item, obs] = bernoulli_logit_rng(mu[item] + lambda[item-1]*theta[obs]);

      // calculate conditional data likelihood for LOO/WAIC
      personLike[obs] = 
        personLike[obs] + 
        bernoulli_logit_lpmf(Y[item, obs] | mu[item] + lambda[item-1]*theta[obs]);
    }
    
  }

}
