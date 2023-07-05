data {
  int<lower=0> nObs;                            // number of observations
  int<lower=0> nItems;                          // number of items
  array[nItems, nObs] int<lower=0, upper=1>  Y; // item responses in an array

  vector[nItems] priorMeanMu;             // prior mean vector for intercept parameters
  matrix[nItems, nItems] priorCovMu;      // prior covariance matrix for intercept parameters

  real priorMeanLambda;
  real priorSDLambda;
}

parameters {
  vector[nObs] theta;                // the latent variables (one for each person)
  vector[nItems] mu;                 // the item intercepts (one for each item)
  real lambda;
}
  
model {
  
  mu ~ multi_normal(priorMeanMu, priorCovMu);   // Prior for item intercepts
  theta ~ normal(0, 1);                         // Prior for latent variable (with mean/sd specified)
  lambda ~ normal(priorMeanLambda, priorSDLambda);                             // Prior for discrimination parameter

  for (item in 1:nItems){
    Y[item] ~ bernoulli_logit(mu[item] + lambda*theta);
  }
  
}

generated quantities{

  // converting item difficulties to item intercepts (easiness)
  vector[nItems] b;
  for (item in 1:nItems){
    b[item] = -1*mu[item]/lambda;
  }
  
  // for PPMC:
  array[nItems, nObs] int<lower=0> simY;
  
  // for LOO/WAIC:
  vector[nObs] personLike = rep_vector(0.0, nObs);
  
  for (item in 1:nItems){

    for (obs in 1:nObs){
      // generate data based on distribution and model
      simY[item, obs] = bernoulli_logit_rng(mu[item] + lambda*theta[obs]);
      
      // calculate conditional data likelihood for LOO/WAIC
      personLike[obs] = 
        personLike[obs] + 
        bernoulli_logit_lpmf(Y[item, obs] | mu[item] + lambda*theta[obs]);
    }
    
  }

}
