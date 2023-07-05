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

  vector[nItems-1] priorMeanLambda;     // prior mean for item discrimination/loadings
  matrix[nItems-1, nItems-1] priorCovLambda;  // prior covariance matrix for item discrimination/loadings

  real priorPersonThetaSDmean;             // prior mean for person theta standard deviation
  real priorPersonThetaSDsd;               // prior standard deviation for person theta standard deviation

  real priorSchoolThetaSDmean;             // prior mean for school theta standard deviation
  real priorSchoolThetaSDsd;               // prior standard deviation for school theta standard deviation
}

parameters {
  array[nSchools] vector[nItems] itemRI;          // the item random intercepts (one for each item for each school)
  vector[nItems] gamma;                           // the item intercepts (one for each item)
  vector<lower=0>[nItems] itemRIsd;                      // standard deviation of item random intercepts
  real<lower=0> personThetaSD;                    // standard deviation of person theta
  real<lower=0> schoolThetaSD;                    // standard deviation of person theta
  vector[nPersons] personTheta;                   // person theta
  vector[nSchools] schoolTheta;                   // school theta
  vector[nItems-1] lambda;                    // item discrimination/loadings
}

model {
  vector[nItems] priorMeanItemRI = rep_vector(0, nItems);
  matrix[nItems, nItems] itemRIcovL;  // lower triangle of cholesky of covariance matrix for item random intercepts

  itemRIsd ~ lognormal(priorItemRIsdMean, priorItemRIsdSD);        // Prior for item random intercept standard deviation
  
  for (item in 1:nItems){
    itemRI[,item] ~ normal(priorMeanItemRI[item], itemRIsd[item]); 
  }
  

  gamma ~ multi_normal(priorMeanGamma, priorCovGamma);  // Prior for item intercepts
  lambda ~ multi_normal(priorMeanLambda, priorCovLambda);  // Prior for item discrimination/loadings

  personThetaSD ~ lognormal(priorPersonThetaSDmean, priorPersonThetaSDsd); // Prior for person theta standard deviation
  personTheta ~ normal(0, personThetaSD);              // Prior for person theta

  schoolThetaSD ~ lognormal(priorSchoolThetaSDmean, priorSchoolThetaSDsd); // Prior for school theta standard deviation
  schoolTheta ~ normal(0, schoolThetaSD);

  for (person in 1:nPersons){
      // first item has withinLambda = 1
      Y[1, person] ~ bernoulli_logit(gamma[1] + personTheta[person] + 
                                     schoolTheta[school[person]] + itemRI[school[person], 1]);
      for (item in 2:(nItems)){
          Y[item, person] ~ bernoulli_logit(gamma[item] + 
                                            lambda[item-1]*(personTheta[person] + schoolTheta[school[person]]) + 
                                            itemRI[school[person], item]);
      }

  }
  
}

generated quantities {

  vector[nPersons] personLike = rep_vector(0.0, nPersons);

  // for PPMC:
  array[nItems, nPersons] int<lower=0> simY;

  for (person in 1:nPersons){
    personLike[person] = 
        personLike[person] + 
        bernoulli_logit_lpmf(Y[1, person] | 
                              gamma[1] + personTheta[person] + schoolTheta[school[person]] + 
                              itemRI[school[person], 1]);
    
    // generate data based on distribution and model
    simY[1, person] = bernoulli_logit_rng(
       gamma[1] + personTheta[person] + schoolTheta[school[person]] + itemRI[school[person], 1]
    );

    for (item in 2:(nItems)){    
      // calculate conditional data likelihood for LOO/WAIC
      personLike[person] = 
        personLike[person] + 
        bernoulli_logit_lpmf(Y[item, person] | 
                                gamma[item] + lambda[item-1]*(personTheta[person] +
                                schoolTheta[school[person]]) + itemRI[school[person], item]);

      // generate data based on distribution and model
      simY[item, person] = bernoulli_logit_rng(
         gamma[item] + lambda[item-1]*(personTheta[person] +
            schoolTheta[school[person]]) + itemRI[school[person], item]
      );
    }
    
  }
}
