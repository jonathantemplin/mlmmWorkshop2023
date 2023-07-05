data {
  int<lower=0> nPersons;                            // number of observations
  int<lower=0> nSchools;                            // number of schools
  int<lower=0> nItems;                          // number of items
  array[nItems, nPersons] int<lower=0, upper=1>  Y; // item responses in an array

  vector[nItems] priorMeanGamma;             // prior mean vector for intercept parameters
  matrix[nItems, nItems] priorCovGamma;      // prior covariance matrix for intercept parameters

  array[nPersons] int school;                   // school number for each person

  vector[nItems] priorItemRIsdMean;         // prior mean for item random intercept standard deviation
  vector[nItems] priorItemRIsdSD;           // prior standard deviation for item random intercept standard deviation

  real priorRIcorrLKJparam;                 // prior parameter for LKJ prior for item random intercept correlation matrix

  vector[nItems] priorMeanWithinLambda;     // prior mean for item discrimination/loadings
  matrix[nItems, nItems] priorCovWithinLambda;  // prior covariance matrix for item discrimination/loadings
}

parameters {
  array[nSchools] vector[nItems] itemRI;          // the item random intercepts (one for each item for each school)
  vector[nItems] gamma;                           // the item intercepts (one for each item)
  cholesky_factor_corr[nItems] itemRIcorrL;       // cholesky factor of the correlation matrix for item random intercepts
  vector<lower=0>[nItems] itemRIsd;                      // standard deviation of item random intercepts

  vector[nPersons] personTheta;                   // person theta
  vector[nItems] withinLambda;                    // item discrimination/loadings
}

model {
  vector[nItems] priorMeanItemRI = rep_vector(0, nItems);
  matrix[nItems, nItems] itemRIcovL;  // lower triangle of cholesky of covariance matrix for item random intercepts

  itemRIsd ~ lognormal(priorItemRIsdMean, priorItemRIsdSD);        // Prior for item random intercept standard deviation
  itemRIcorrL ~ lkj_corr_cholesky(priorRIcorrLKJparam);            // Prior for item random intercept correlation matrix
  itemRIcovL = diag_pre_multiply(itemRIsd, itemRIcorrL);           // Form lower triangle of Covariance matrix for item random intercepts
  itemRI ~ multi_normal_cholesky(priorMeanItemRI, itemRIcovL); 

  gamma ~ multi_normal(priorMeanGamma, priorCovGamma);  // Prior for item intercepts
  withinLambda ~ multi_normal(priorMeanWithinLambda, priorCovWithinLambda);  // Prior for item discrimination/loadings

  personTheta ~ normal(0, 1);              // Prior for person theta

  for (person in 1:nPersons){
      for (item in 1:nItems){
          Y[item, person] ~ bernoulli_logit(gamma[item] + withinLambda[item]*personTheta[person] + itemRI[school[person], item]);
      }
  }
  
}

generated quantities {
  corr_matrix[nItems] itemRIcorr;  // correlation matrix for item random intercepts
  cholesky_factor_cov[nItems] itemRIcovL;  // cholesky factor of covariance matrix for item random intercepts
  cov_matrix[nItems] itemRIcov;  // covariance matrix for item random intercepts

  itemRIcorr = multiply_lower_tri_self_transpose(itemRIcorrL);  // correlation matrix for item random intercepts
  itemRIcovL = diag_pre_multiply(itemRIsd, itemRIcorrL);  // cholesky factor of covariance matrix for item random intercepts
  itemRIcov = multiply_lower_tri_self_transpose(itemRIcovL);  // covariance matrix for item random intercepts

  vector[nPersons] personLike = rep_vector(0.0, nPersons);
  
  // for PPMC:
  array[nItems, nPersons] int<lower=0> simY;

  for (person in 1:nPersons){
    for (item in 1:nItems){    
      // calculate conditional data likelihood for LOO/WAIC
      personLike[person] = 
        personLike[person] + 
        bernoulli_logit_lpmf(Y[item, person] | gamma[item] + 
          withinLambda[item]*personTheta[person] + itemRI[school[person], item]);
      
      // generate data based on distribution and model
      simY[item, person] = bernoulli_logit_rng(
        gamma[item] + withinLambda[item]*personTheta[person] + itemRI[school[person], item]);
    }
  }
}

