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

  real priorSlopePersonThetaFRlunchMean;  // prior mean for slope of person theta on frlunch
  real priorSlopePersonThetaFRlunchSD;    // prior standard deviation for slope of person theta on frlunch

  vector[2] priorSchoolRIsdMean;
  vector[2] priorSchoolRIsdSD;

  real priorSchoolRIcorrLKJparam;               // prior for item random intercept correlation matrix
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

  cholesky_factor_corr[2] schoolRIcorrL;       // cholesky factor of the correlation matrix for item random intercepts
  vector<lower=0>[2] schoolRIsd;                      // standard deviation of item random intercepts

  // parameters for frlunch model
  real frlunchIntercept;
  array[nSchools] vector[2] schoolRI;
  real<lower=0> frlunchRIsd;
  real slopePersonTheta_Frlunch;
}

model {

  vector[2] priorMeanSchoolRI = rep_vector(0, 2);
  matrix[2, 2] schoolRIcovL;  // lower triangle of cholesky of covariance matrix for item random intercepts

  // random intercepts for school (both frlunch and item random intercept factor)
  schoolRIsd ~ lognormal(priorSchoolRIsdMean, priorSchoolRIsdSD);        // Prior for item random intercept standard deviation
  schoolRIcorrL ~ lkj_corr_cholesky(priorSchoolRIcorrLKJparam);            // Prior for item random intercept correlation matrix
  schoolRIcovL = diag_pre_multiply(schoolRIsd, schoolRIcorrL);           // Form lower triangle of Covariance matrix for item random intercepts
  schoolRI ~ multi_normal_cholesky(priorMeanSchoolRI, schoolRIcovL); 

  // priors for predicting thetas at person level
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

  for (person in 1:nPersons){
    personTheta[person] ~ normal(0 + slopePersonTheta_Frlunch*frlunch[person], personThetaSD);   
    frlunch[person] ~ bernoulli_logit(frlunchIntercept + schoolRI[school[person],1]);
    
    // first item has withinLambda = 1
    Y[1, person] ~ bernoulli_logit(gamma[1] + 
                                   personTheta[person] + schoolRI[school[person],2] + 
                                   itemRI[school[person], 1]);
    


    for (item in 2:(nItems)){
        Y[item, person] ~ bernoulli_logit(gamma[item] + 
                                          lambda[item-1]*(personTheta[person] + schoolRI[school[person],2]) + 
                                          itemRI[school[person], item]);
    }

  }
  
}

generated quantities {

  corr_matrix[2] schoolRIcorr;  // correlation matrix for item random intercepts
  cholesky_factor_cov[2] schoolRIcovL;  // cholesky factor of covariance matrix for item random intercepts
  cov_matrix[2] schoolRIcov;  // covariance matrix for item random intercepts

  schoolRIcorr = multiply_lower_tri_self_transpose(schoolRIcorrL);  // correlation matrix for item random intercepts
  schoolRIcovL = diag_pre_multiply(schoolRIsd, schoolRIcorrL);  // cholesky factor of covariance matrix for item random intercepts
  schoolRIcov = multiply_lower_tri_self_transpose(schoolRIcovL);  // covariance matrix for item random intercepts
   
  real schoolThetaSlope;

  schoolThetaSlope = schoolRIcov[1,2]/schoolRIcov[1,1];

  real residualSchoolThetaVar;
  real schoolPseudoR2;

  residualSchoolThetaVar = schoolRIcov[2,2] - square(schoolRIcov[1,2])/schoolRIcov[1,1];
  schoolPseudoR2 = (schoolRIcov[2,2] - residualSchoolThetaVar)/schoolRIcov[2,2];
}
