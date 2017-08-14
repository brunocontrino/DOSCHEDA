#' An S4 class to run the doscheda pipeline
#'
#' @slot input A data.frame containing the input data
#' @slot normData A data.frame containin a processed and standardised version of the input data
#' @slot finalData A data.frame containing the final data produced by the pipline
#' @slot parameters A list containing all the parameters required to make the pipeline run successfully
#' @slot datasets A list containing other potentially useful datasets
#'
#' @importFrom methods new
#' @exportClass ChemoProtSet
#' @export
#'
DoschedaSet <- setClass("ChemoProtSet", slots = c(input = "data.frame", normData = "data.frame", 
    finalData = "data.frame", parameters = "list", datasets = "list"), prototype = list(input = data.frame(), 
    normData = data.frame(), finalData = data.frame(), parameters = list(chans = 3, reps = 2, chanNames = NA, 
        PD = TRUE, dataType = "LFC", modelType = "linear", removePeps = NA, organism = "h.sapiens", 
        pearsonThrsh = 0.4), datasets = list(GeneID = data.frame())), validity = function(object) {
    if (!is.data.frame(object@input)) {
        return("Input data set is not a data.frame, please input a data.frame")
    } else if (!is.data.frame(object@normData)) {
        return("Normalised data set is not a data.frame, please input a data.frame")
    } else if (!is.data.frame(object@finalData)) {
        return("final data set is not a data.frame, please input a data.frame")
    }
    return(TRUE)
})



