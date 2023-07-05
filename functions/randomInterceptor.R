#' Estimate empty models with and without random intercepts
#'
#' @description
#' Estimates empty models with and without random intercepts across packages and reports ICCs.
#' Calculates intraclass correlation and reports likelihood ratio test statistic, degrees of
#' freedom, and p-value.
#'
#' @param data A data frame containing numeric variables and a unit ID variable
#' @param outcome A character variable containing the name of the outcome variable from \code{data}
#' @param level2ID A character variable containing the name of the level-two ID variable from \code{data}
#' @param REML A logical value where \code{TRUE} (default) estimates the models using REML and
#'             \code{FALSE} estimates the models using ML
#'
#' @return A list with all calculated values along with a display of these values to the console window
#'
#' @examples
#'
#' # load data for example 3a
#' data(example3aData)
#'
#' # remove all NA values from langpost
#' newData = example3aData[complete.cases(example3aData$langpost),]
#'
#' # estimate models using REML with langpost DV and schoolID as level-two ID variable:
#' initLangpost = randomInterceptor(data = example3aData,
#'                                  outcome = "langpost",
#'                                  level2ID = "schoolID",
#'                                  REML = TRUE)
#' # OUTPUT
#' # Empty Level-1 vs Empty Level-2 Random Intercept Model Comparison
#' # ICC: 0.2191
#' # Test Statistic: 502.90794
#' # Test df: 1
#' # Test p-value: 0
#'
#' @export
randomInterceptor = function(data, outcome, level2ID, REML=TRUE){
  # estimate gls() for model w/o random intercept
  if (!require(nlme) | !require(lme4)) stop("function needs nlme and lme4 packages installed; please install and rerun")

  if (REML){
    methodName = "REML"
  } else{
    methodName = "ML"
  }
  baseModel = formula(paste(outcome, "~", 1))
  glsModel = nlme::gls(model = baseModel, data = data, method = methodName)
  summaryGlsModel = summary(glsModel)

  RImodel = formula(paste0(outcome, "~ 1 + (1|", level2ID, ")"))
  lmerModel = lme4::lmer(formula = RImodel, data = data, REML = REML)
  summaryLmerModel = summary(lmerModel)

  tao0_2 = summaryLmerModel$varcor[[1]][1]
  sigma_2 = summaryLmerModel$sigma^2
  ICC = tao0_2/(tao0_2 + sigma_2)

  LL_LmerModel = summaryLmerModel$logLik[1]
  LL_glsModel = summaryGlsModel$logLik

  LRTstatistic = -2*(LL_glsModel - LL_LmerModel)
  LRTdf = 1
  LRTpvalue = 1-pchisq(q = LRTstatistic, df = 1)

  cat(paste("\nEmpty Level-1 vs Empty Level-2 Random Intercept Model Comparison\n\n"))
  cat(paste("ICC:", round(ICC,5), "\n"))
  cat(paste("Test Statistic:", round(LRTstatistic, 5), "\n"))
  cat(paste("Test df:", LRTdf, "\n"))
  cat(paste("Test p-value:", round(LRTpvalue, 5), "\n"))

  return(list(ICC = ICC, modelLL_level1 = LL_glsModel, modelLL_RI = LL_LmerModel,
              LRTstatistic = LRTstatistic, LRTpvalue  = LRTpvalue, LRTdf = LRTdf))
}
