#' Provides marginal R squared for comparing a model with an empty model
#'
#' @description
#' Provides marginal R squared for comparing a model with an empty model
#'
#' @param model A model (from \code{lme4::lmer}) that all predictors from \code{largerModel}
#' @param dvName A character variable containing the name of the outcome variable from \code{data}
#' @param data A data frame containing numeric variables and a unit ID variable
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
#' # add school-level homework variable as predictor
#' hwModel = lmer(formula = langpost ~ 1 + homework + (1|schoolID), data = newData, REML = TRUE)
#' summary(hwModel)
#'
#' # calculate total (marginal) R squared using this function
#' hwVSemptyTotalR2 = totalRSquaredinator(model = hwModel, dvName = "langpost", data = newData)
#'
#' # OUTPUT
#' # Total R2 Estimate
#' # Total R2: 0.00126
#'
#' @export
totalRSquaredinator = function(model, dvName, data){
  totalVariance = as.numeric(cor(data[dvName], predict(model, re.form=NA))^2)
  names(totalVariance) = "total R2"
  cat("\nTotal R2 Estimate\n")
  cat(paste0("Total R2: ", round(totalVariance[1],5)))
  return(totalVariance)
}
