# Multilevel Modeling: Predicting Observed Sum Score for Students Nested in Schools 
# Using ML Estimation Below

# Preliminary Steps

# Set width of output and number of significant digits printed,
# number of digits before using scientific notation, shut off significance stars, # cores
options(width=120, digits=8, scipen=9, show.signif.stars=FALSE, mc.cores=4)

#####  Check to see if packages are downloaded, install if not, then load  #####

# To get compact data description
if (!require("psych")) install.packages("psych"); library(psych) 

# To estimate MLMs using gls or lme
if (!require("nlme")) install.packages("nlme"); library(nlme) 

# To estimate MLMs using lmer
if (!require("lme4")) install.packages("lme4"); library(lme4) 

# To get Satterthwaite DDF in lmer
if (!require("lmerTest")) install.packages("lmerTest"); library(lmerTest) 

# To get ICC in lmer
if (!require("performance")) install.packages("performance"); library(performance)

# To estimate multivariate MLM using multilevel SEM
if (!require("lavaan")) install.packages("lavaan"); library(lavaan)

# Clear workspace (re-run as needed for troubleshooting purposes)
rm(list = ls())

###############################################################################################
# Data Import, Manipulation, and Description

# Set working directory
setwd("C:/Dropbox/mlmmWorkshop2023")

# Load R functions for this example from folder within working directory
functions = paste0("functions/",dir("functions/"))
temp = lapply(X = functions, FUN = source)

# Load example data from folder within working directory
load("data/modelingData.RData")

# Create school means for sum score student outcome and free/reduced lunch student predictor 
modelingData = addUnitMeans(data=modelingData, unitVariable="schoolID", 
                            meanVariables=c("sumScore","frlunch"), 
                            newNames=c("SMsumScore","SMfrlunch"))

# Descriptive statistics for student variables and new school means
print(describe(modelingData[c("sumScore","frlunch","SMsumScore","SMfrlunch")]),digits=3)

# Constant-center school lunch near mean to use as observed level-2 predictor
modelingData$SMfrlunch30 = modelingData$SMfrlunch - .30

# Cluster-mean-center student lunch at school mean to use as observed level-1 predictor
modelingData$WSfrlunch = modelingData$frlunch - modelingData$SMfrlunch

#################################################################################################
# Models for Partitioning Student-Level from School-Level Variance

# Partitioning Variance in the Sum Score Outcome using General Models

# Single-level empty model predicting observed sum score ignoring school
# Using gls instead of lm to get model log-likelihood for model comparison
modelEmptyGLM = gls(data=modelingData, method="ML", model=sumScore~1)
summary(modelEmptyGLM)

# Two-level empty model predicting observed sum score with students nested in schools
modelEmptyRI = lmer(data=modelingData, REML=FALSE, sumScore~1+(1|schoolID))
summary(modelEmptyRI)

# Show intraclass correlation and its likelihood ratio test
icc(modelEmptyRI); ranova(modelEmptyRI)


# Partitioning Variance in the Binary Free/Reduced Lunch Predictor using Generalized Models

# Single-level empty model predicting observed free/reduced lunch ignoring school
modelEmptyGLMfr = glm(data=modelingData, family=binomial(link="logit"), formula=frlunch~1)
summary(modelEmptyGLMfr) # Null deviance= -2LL already
# Convert logit intercept into probability
modelEmptyGLMfrProb=1/(1+exp(-1*coefficients(modelEmptyGLMfr))); modelEmptyGLMfrProb 

# Two-level empty model predicting observed free/reduced lunch with students nested in schools
modelEmptyRIfr = glmer(data=modelingData, family=binomial(link="logit"), frlunch~1+(1|schoolID))
summary(modelEmptyRIfr) # deviance = -2LL already

# Convert logit intercept into probability (sub-object beta holds fixed intercept)
modelEmptyRIfrProb=1/(1+exp(-1*(modelEmptyRIfr@beta))); modelEmptyRIfrProb 

# Compute ICC using pi^2/3 = 3.29 as residual variance (sub-object theta holds random intercept variance)
modelEmptyRIfr@theta^2/(modelEmptyRIfr@theta^2+(pi^2/3)) 

# Likelihood Ratio Test for Addition of Random Intercept Variance
DevTest=-2*(logLik(modelEmptyGLMfr)-logLik(modelEmptyRIfr))
Pvalue=pchisq((DevTest), df=1, lower.tail=FALSE)
# Test Statistic and P-values for DF=1 
DevTest; Pvalue


##################################################################################################
#  Models Predicting the Observed Sum Score from Free/Reduced Lunch for Students Nested in Schools

# Add smushed level-1 slope for frlunch
modelSmushed = lmer(data=modelingData, REML=FALSE, sumScore~1+frlunch+(1|schoolID))
summary(modelSmushed, ddf="Satterthwaite")

# Proportion explained of each variance component relative to empty model
pseudoRSquaredinator(smallerModel=modelEmptyRI, largerModel=modelSmushed)

# Add centered school mean to unsmush the level-1 slope for frlunch
modelLunch = lmer(data=modelingData, REML=FALSE, sumScore~1+frlunch+SMfrlunch30+(1|schoolID))
summary(modelLunch, ddf="Satterthwaite")

# Compute full between level-2 effect
contest1D(modelLunch, L=c(0,1,1))

# Proportion explained of each variance component relative to smushed model
pseudoRSquaredinator(smallerModel=modelSmushed, largerModel=modelLunch)

# Proportion explained of each variance component relative to empty model
pseudoRSquaredinator(smallerModel=modelEmptyRI, largerModel=modelLunch)

# Total R2 relative to empty model
totalRSquaredinator(model=modelLunch, dvName="sumScore", data=modelingData)

# Cluster-mean-centered version of modelLunch
modelLunchCMC = lmer(data=modelingData, REML=FALSE, sumScore~1+WSfrlunch+SMfrlunch30+(1|schoolID))
summary(modelLunchCMC, ddf="Satterthwaite")

# Compute contextual level-2 effect
contest1D(modelLunchCMC, L=c(0,-1,1))

# Total R2 relative to empty model
totalRSquaredinator(model=modelLunchCMC, dvName="sumScore", data=modelingData)


##################################################################################################
# Multivariate MLM using Latent Centering, Pretending frlunch is Continuous

MultivSyntax = "
level: 1
  # Level-1 residual variance for SumScore only
    sumScore ~~ sumScore
  # frlunch predicts sumScore: level-1 within-school slope
    sumScore ~ (within)*frlunch
level: 2
  # Fixed intercepts
    sumScore ~ 1; frlunch ~ 1
  # Level-2 random intercept variances
    sumScore ~~ sumScore
    frlunch ~~ frlunch
  # frlunch predicts sumScore: level-2 contextual slope
    sumScore ~ (context)*frlunch
  # Compute full between level-2 effect
    between := within + context
"
modelMultiv = lavaan(model=MultivSyntax, data=modelingData, cluster="schoolID", 
                     mimic="mplus", std.lv=FALSE, estimator="ML") 
summary(object=modelMultiv)


##################################################################################################
#  Models with a Random Slope across Schools for Cluster-Mean-Centered Student Free/Reduced Lunch 

# Add random level-1 slope for CMC frlunch
modelRandSlope = lmer(data=modelingData, REML=FALSE, 
                      sumScore~1+WSfrlunch+SMfrlunch30+(1+WSfrlunch|schoolID))
summary(modelRandSlope, ddf="Satterthwaite")
# Likelihood ratio test for significance of random slope variance
ranova(modelRandSlope)

# Add cross-level interaction predicting random level-1 slope for CMC frlunch (demo purposes only)
modelCrossLevel = lmer(data=modelingData, REML=FALSE, 
                       sumScore~1+WSfrlunch+SMfrlunch30+WSfrlunch:SMfrlunch30+(1+WSfrlunch|schoolID))
summary(modelCrossLevel, ddf="Satterthwaite")


##################################################################################################
# Example of how to export a .csv file for use in Mplus
# Copy data, replace all missing values with -999 for Mplus
modelingData_Mplus = modelingData
modelingData[is.na(modelingData)] = -999

# Write to .csv file with column names (to be removed later)
write.table(x=modelingData_Mplus, col.names=TRUE, row.names=FALSE, sep=",",
            file="modelingData_names.csv")

