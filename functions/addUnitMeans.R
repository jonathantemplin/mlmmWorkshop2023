#' Add cluster/group/person means to a data frame for a set of specified variables
#'
#' @description
#' Add cluster/group/person means to a data frame
#'
#' @param data A data frame containing numeric variables and a unit ID variable
#' @param unitVariable A variable containing the ID number for each cluster/group/person
#' @param meanVariables A character vector containing the names of the variables for which
#'                      the means for each cluster/group/person are needed
#' @param newNames A character vector of the same length as \code{meanVariables} with the names
#'                 for the new variables of cluster/group/person means
#'
#' @return A data frame where new variables are added
#'
#' @examples
#'
#' # load data for example 3a
#' data(example3aData)
#'
#' # display variable names for data set
#' names(example3aData)
#'
#' # add school means with this function
#' analysisData = addUnitMeans(
#'   data=studentsZ,
#'   unitVariable = "schoolID",
#'   meanVariables = c("langpost", "IQverbz", "IQperfz"),
#'   newNames = c("SMlangpost", "SMIQverbz", "SMIQperfz")
#' )
#'
#' # display new names
#' names(analysisData)
#'
#' @export
addUnitMeans = function(data, unitVariable, meanVariables, newNames){
  # create unit-level data set:

  unitMeans = t(sapply(
    X = unique(data[,unitVariable]),
    FUN = function(x, data) {
      return(c(
        x,
        length(which(data[,unitVariable] == x)),
        apply(
          X = as.data.frame(data[which(data[,unitVariable] == x), meanVariables]), 
          MARGIN = 2, 
          FUN = mean, 
          rm.na=TRUE
          )
        )
      )
    },data = data
  ))

  unitMeans = data.frame(unitMeans)
  names(unitMeans) = c(unitVariable, paste0("Nper", unitVariable), newNames)

  newData = merge(x = data, y = unitMeans, by = unitVariable)
  return(newData)
}


