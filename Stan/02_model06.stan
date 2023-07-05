data {
  int<lower=0> N;                            // total number of observations
  int<lower=0> nSchools;                     // number of unique level-2 units in data
  array[N] int obsSchoolID;                  // the level-2 ID corresponding to each row of the data
  array[N] int frlunch;                      // person-level frlunch variable
  vector[N] score;                           // score outcome variable
  real priorBetaMean;                        // prior mean vector for coefficients (all same)
  real priorBetaSD;                          // prior sd for coefficients (all same)
  real priorSigmaMean;                       // prior mean parameter for residual standard deviation
  real priorSigmaSD;                         // prior sd parameter for residual standard deviation
  real priorTauScoreMean;                    // prior mean parameter for score random intercept standard deviation
  real priorTauScoreSD;                      // prior sd parameter for score random intercept standard deviation
  real priorTauFrlunchMean;                  // prior mean parameter for frlunch random intercept standard deviation
  real priorTauFrlunchSD;                    // prior sd parameter for frlunch random intercept standard deviation
}

parameters {
  real scoreIntercept;                       // fixed intercept for score
  real frlunchIntercept;                     // fixed intercept for frlunch
  real frlunchSlope;                         // fixed slope for frlunch
  real frlunchMeanSlope;                     // fixed slope for frlunch mean (from Multivariate MLM)
  real<lower=0> sigma;                       // residual standard deviation for score
  vector[nSchools] scoreRandomIntercept;     // random intercept for each level-2 unit for score
  vector[nSchools] frlunchRandomIntercept;   // random intercept for each level-2 unit for frlunch
  real<lower=0> tauScoreIntercept;           // random intercept standard deviation for score
  real<lower=0> tauFrlunchIntercept;         // random intercept standard deviation for frlunch
}

model {
  scoreIntercept ~ normal(priorBetaMean, priorBetaSD);                     // prior for score intercept
  frlunchIntercept ~ normal(priorBetaMean, priorBetaSD);                   // prior for frlunch intercept
  frlunchSlope ~ normal(priorBetaMean, priorBetaSD);                       // prior for frlunch slope
  frlunchMeanSlope ~ normal(priorBetaMean, priorBetaSD);                   // prior for frlunch mean slope

  sigma ~ lognormal(priorSigmaMean, priorSigmaSD);                         // prior for sigma
  tauScoreIntercept ~ lognormal(priorTauScoreMean, priorTauScoreSD);       // prior for tau for score
  tauFrlunchIntercept ~ lognormal(priorTauFrlunchMean, priorTauFrlunchSD); // prior for tau for frlunch

  for (school in 1:nSchools){
    scoreRandomIntercept[school] ~ normal(0, tauScoreIntercept);
    frlunchRandomIntercept[school] ~ normal(0, tauFrlunchIntercept);
  }

  for (obs in 1:N){
    frlunch[obs] ~ bernoulli_logit(frlunchIntercept + frlunchRandomIntercept[obsSchoolID[obs]]);

    score[obs] ~ normal(scoreIntercept + frlunchSlope*frlunch[obs] + 
                    frlunchMeanSlope * (frlunchIntercept + frlunchRandomIntercept[obsSchoolID[obs]]) + 
                    scoreRandomIntercept[obsSchoolID[obs]], sigma); // linear model
    
  }
  
}

generated quantities {
  real ICCscore;
  ICCscore = tauScoreIntercept^2/(tauScoreIntercept^2 + sigma^2);

  real ICCfrlunch;
  ICCfrlunch = tauFrlunchIntercept^2/(tauFrlunchIntercept^2 + (pi()^2/3));

  vector[N] personLike = rep_vector(0.0, N);
  for (obs in 1:N){
    personLike[obs] = 
      normal_lpdf(score[obs] | scoreIntercept + frlunchSlope*frlunch[obs] + 
                  frlunchMeanSlope * (frlunchIntercept + frlunchRandomIntercept[obsSchoolID[obs]]) + 
                  scoreRandomIntercept[obsSchoolID[obs]], sigma);
  }
}
