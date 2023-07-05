#' Provides pseudo R squared for comparing two models (one nested within the other)
#'
#' @description
#' Provides pseudo R squared for comparing two models (one nested within the other).
#'
#' @param smallerModel A model (from \code{lme4::lmer}) that has a subset of predictors from \code{largerModel}
#' @param largerModel A  model (from \code{lme4::lmer}) that has all predictors from \code{smallerModel} and additional predictors
#'
#' @return A list with all calculated values along with a display of these values to the console window
#'
#' @examples
#'
#' # load data for example 3a
#' data(example3aData)
#'
#' # remove all NA values from langpost and
#' newData = example3aData[complete.cases(example3aData[,c("langpost", "homework")]),]
#'
#' # estimate empty model with random intercept for school
#' emptyModel = lmer(formula = langpost ~ 1 + (1|schoolID), data = newData, REML = TRUE)
#' summary(emptyModel)
#'
#' # add school-level homework variable as predictor
#' hwModel = lmer(formula = langpost ~ 1 + homework + (1|schoolID), data = newData, REML = TRUE)
#' summary(hwModel)
#'
#' # calculate pseudo R squared using this function
#' hwVSempty = pseudoRSquaredinator(smallerModel = emptyModel, largerModel = hwModel)
#'
#' # OUTPUT
#' # Pseudo R2 Estimates
#' # R2 Random.(Intercept): 0.00179
#' # R2 L1.sigma2: 5e-05
#'
#' @export
pseudoRSquaredinator = function(smallerModel, largerModel){
  if (!require(lme4)) stop("Please install lme4 before using this function.")
  if (length(grep(x = class(smallerModel), pattern="lmerMod"))==0) stop("Models must be from lmer. Problem with smallerModel.")
  if (length(grep(x = class(largerModel), pattern="lmerMod"))==0) stop("Models must be from lmer. Problem with largerModel.")

  # get model summaries
  smallModelSummary = summary(smallerModel)
  largeModelSummary = summary(largerModel)

  # number of variances
  nRandomEffectsSmall = nrow(smallModelSummary$varcor[[1]])
  nRandomEffectsLarge = nrow(largeModelSummary$varcor[[1]])

  if (nRandomEffectsSmall != nRandomEffectsLarge) stop("Models have different numbers of random effects. Please ensure both models have same random effects.")

  smallG = smallModelSummary$varcor[[1]][1:nRandomEffectsSmall, 1:nRandomEffectsSmall]
  largeG = largeModelSummary$varcor[[1]][1:nRandomEffectsLarge, 1:nRandomEffectsLarge]

  pseudoR2 = NULL
  namesVec = NULL

  if (nRandomEffectsLarge == 1){ # single random effect appears as scalar not matrix
    pseudoR2 = c(pseudoR2, (smallG[1]-largeG[1])/smallG[1])
    namesVec = c(namesVec, paste0("R2 Random.",colnames(smallModelSummary$varcor[[1]])[1]))
  } else {
    for (re in 1:nRandomEffectsLarge){
      pseudoR2 = c(pseudoR2, (smallG[re,re]-largeG[re,re])/smallG[re,re])
      namesVec = c(namesVec, paste0("R2 Random.",colnames(smallModelSummary$varcor[[1]])[re]))
    }
  }


  # get sigma2
  sigma2small = smallModelSummary$sigma^2
  sigma2large = largeModelSummary$sigma^2
  pseudoR2 = c(pseudoR2, (sigma2small-sigma2large)/sigma2small)

  namesVec = c(namesVec, "R2 L1.sigma2")
  names(pseudoR2) = namesVec

  cat("\nPseudo R2 Estimates\n")
  for (p in 1:length(pseudoR2)){
    cat(paste0(names(pseudoR2)[p], ": ", round(pseudoR2[p], 5), "\n"))
  }

  return(pseudoR2)
}
