#' MeanSd plot for objects of class ChemoProtSet
#'
#' Shows the ranked means with a running median calculated with a window size of 10%
#'
#' @param x object of class 'ChemoProtSet'
#' @param ... other plot options
#'
#' @return meanSd plot for objects of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname meanSdPlot-methods
#' @examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'meanSdPlot(ex)

setGeneric(name="meanSdPlot",
           def=function(x, ...)
           {
             standardGeneric("meanSdPlot")
           }
)

#' @rdname meanSdPlot-methods
#' @importFrom vsn meanSdPlot
#' @import grDevices
#' @aliases meanSdPlot,ANY,ANY-method
#'
setMethod(f= "meanSdPlot",
          signature="ChemoProtSet",
          definition=function(x, ...)
          {
            if(x@parameters$dataType =='linear'){
              vsn::meanSdPlot(as.matrix(x@normData[,1:(x@parameters$chans * x@parameters$reps)]))
            } else{
              vsn::meanSdPlot(as.matrix(x@normData[,1:((x@parameters$chans - 1)* x@parameters$reps)]))
            }
          }

)
