data {
  int<lower=0> nObs;                            // number of observations
  int<lower=0> nItems;                          // number of items
  array[nItems, nObs] int<lower=0, upper=1>  Y; // item responses in an array

  vector[nItems] priorMeanMu;             // prior mean vector for intercept parameters
  matrix[nItems, nItems] priorCovMu;      // prior covariance matrix for intercept parameters

  vector[nItems] priorMeanLambda;         // prior mean vector for loading parameters
  matrix[nItems, nItems] priorCovLambda;   // prior covariance matrix for loading parameters
}

parameters {
  vector[nObs] theta;                // the latent variables (one for each person)
  vector[nItems] mu;                 // the item intercepts (one for each item)
  vector[nItems] lambda;             // the factor loadings/item discriminations (one for each item)
}

model {

  lambda ~ multi_normal(priorMeanLambda, priorCovLambda); // Prior for item discrimination/factor loadings
  mu ~ multi_normal(priorMeanMu, priorCovMu);             // Prior for item intercepts
  
  theta ~ normal(0, 1);                         // Prior for latent variable (with mean/sd specified)
  
  for (item in 1:nItems){
    Y[item] ~ bernoulli_logit(mu[item] + lambda[item]*theta);
  }
  
}


generated quantities{

  // converting item difficulties to item intercepts (easiness)
  vector[nItems] b;

  for (item in 1:nItems){
    b[item] = -1*mu[item]/lambda[item];
  }
  
  // for PPMC:
  array[nItems, nObs] int<lower=0> simY;
  
  // for LOO/WAIC:
  vector[nObs] personLike = rep_vector(0.0, nObs);
  
  for (item in 1:nItems){

    for (obs in 1:nObs){
      // generate data based on distribution and model
      simY[item, obs] = bernoulli_logit_rng(mu[item] + lambda[item]*theta[obs]);

      // calculate conditional data likelihood for LOO/WAIC
      personLike[obs] = 
        personLike[obs] + 
        bernoulli_logit_lpmf(Y[item, obs] | mu[item] + lambda[item]*theta[obs]);
    }
    
  }

}
