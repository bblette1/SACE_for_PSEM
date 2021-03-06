#------------------------------------------------------------------------------#
#' Check an array object
#'
#' Handles the case where there is a single estimating equation. This function
#' assumes that the object
#'
#' @param object the object to check whether it is an array
#' @return an array - either the orginal object or the given object converted
#' to an array
#------------------------------------------------------------------------------#

check_array <- function(object){
  if(!is.array(object)){
    stopifnot(is.numeric(object))
    array(object, dim = c(1, 1, length(object)))
  } else {
    object
  }
}