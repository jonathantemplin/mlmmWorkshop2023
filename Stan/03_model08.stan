data {
  int<lower=0> nPersons;                            // number of observations
  int<lower=0> nSchools;                            // number of schools
  int<lower=0> nItems;                          // number of items
  array[nItems, nPersons] int<lower=0, upper=1>  Y; // item responses in an array

  vector[nItems] priorMeanGamma;             // prior mean vector for intercept parameters
  matrix[nItems, nItems] priorCovGamma;      // prior covariance matrix for intercept parameters

  array[nPersons] int school;                   // school number for each person
  array[nPersons] int frlunch;                  // frlunch number for each person

  vector[nItems] priorItemRIsdMean;         // prior mean for item random intercept standard deviation
  vector[nItems] priorItemRIsdSD;           // prior standard deviation for item random intercept standard deviation

  vector[nItems-1] priorMeanLambda;     // prior mean for item discrimination/loadings
  matrix[nItems-1, nItems-1] priorCovLambda;  // prior covariance matrix for item discrimination/loadings

  real priorPersonThetaSDmean;             // prior mean for person theta standard deviation
  real priorPersonThetaSDsd;               // prior standard deviation for person theta standard deviation

  real priorSchoolThetaSDmean;             // prior mean for school theta standard deviation
  real priorSchoolThetaSDsd;               // prior standard deviation for school theta standard deviation

  real priorFrlunchRIsdMean;              // prior mean for frlunch random intercept standard deviation
  real priorFrlunchRIsdSD;                // prior standard deviation for frlunch random intercept standard deviation

  real priorSlopeSchoolThetaFRlunchMean;  // prior mean for slope of school theta on frlunch
  real priorSlopeSchoolThetaFRlunchSD;    // prior standard deviation for slope of school theta on frlunch

  real priorSlopePersonThetaFRlunchMean;  // prior mean for slope of person theta on frlunch
  real priorSlopePersonThetaFRlunchSD;    // prior standard deviation for slope of person theta on frlunch

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

  // parameters for frlunch model
  real frlunchIntercept;
  vector[nSchools] frlunchRandomIntercept;
  real<lower=0> frlunchRIsd;
  real slopeSchoolTheta_Frlunch;
  real slopePersonTheta_Frlunch;
}

model {
  // Prior for frlunch random intercept standard deviation
  frlunchRIsd ~ lognormal(priorFrlunchRIsdMean, priorFrlunchRIsdSD);       
  frlunchRandomIntercept ~ normal(0, frlunchRIsd); // Prior for frlunch random intercepts

  // priors for predicting thetas at person and school level
  slopeSchoolTheta_Frlunch ~ normal(priorSlopeSchoolThetaFRlunchMean, priorSlopeSchoolThetaFRlunchSD);
  slopePersonTheta_Frlunch ~ normal(priorSlopePersonThetaFRlunchMean, priorSlopePersonThetaFRlunchSD);

  vector[nItems] priorMeanItemRI = rep_vector(0, nItems);
  matrix[nItems, nItems] itemRIcovL;  // lower triangle of cholesky of covariance matrix for item random intercepts

  itemRIsd ~ lognormal(priorItemRIsdMean, priorItemRIsdSD);        // Prior for item random intercept standard deviation
  
  for (item in 1:nItems){
    itemRI[,item] ~ normal(priorMeanItemRI[item], itemRIsd[item]); 
  }
  

  gamma ~ multi_normal(priorMeanGamma, priorCovGamma);  // Prior for item intercepts
  lambda ~ multi_normal(priorMeanLambda, priorCovLambda);  // Prior for item discrimination/loadings

  personThetaSD ~ lognormal(priorPersonThetaSDmean, priorPersonThetaSDsd); // Prior for person theta standard deviation
  schoolThetaSD ~ lognormal(priorSchoolThetaSDmean, priorSchoolThetaSDsd); // Prior for school theta standard deviation

  for (schoolID in 1:nSchools){
    schoolTheta[school] ~ normal(0 +
     slopeSchoolTheta_Frlunch * (frlunchRandomIntercept[schoolID]), schoolThetaSD); 
  }

  for (person in 1:nPersons){
    personTheta[person] ~ normal(0 + slopePersonTheta_Frlunch*frlunch[person], personThetaSD);   
    frlunch[person] ~ bernoulli_logit(frlunchIntercept + frlunchRandomIntercept[school[person]]);
    
    // first item has withinLambda = 1
    Y[1, person] ~ bernoulli_logit(gamma[1] + 
                                   personTheta[person] + schoolTheta[school[person]] + 
                                   itemRI[school[person], 1]);
    


    for (item in 2:(nItems)){
        Y[item, person] ~ bernoulli_logit(gamma[item] + 
                                          lambda[item-1]*(personTheta[person] + schoolTheta[school[person]]) + 
                                          itemRI[school[person], item]);
    }

  }
  
}

