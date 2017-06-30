#' Method to remove peptides from input data of an object of class 'ChemoProtSet'
#'
#' @param x object of class 'ChemoProtSet'
#' @param normalise string indicating the type of normalisation that should take place ('loess', 'median', 'none')
#'
#' @return  object of class ChemoProtSet
#'
#'@examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex
#'
#' @export
#' @docType methods
#' @rdname runNormalisation-methods
setGeneric(name="runNormalisation",
           def=function(x, normalise = 'loess')
           {
             standardGeneric("runNormalisation")
           }
)
#' @rdname runNormalisation-methods
#' @aliases runNormalisation,ANY,ANY-method
#'
setMethod(f="runNormalisation",
          signature="ChemoProtSet",
          definition=function(x, normalise = 'loess')
          {
             if(x@parameters$dataType  == 'intensity'){
              x@normData <- normalize_data( dataFrame = x@normData,chans = x@parameters$chans,
                                                    reps = x@parameters$reps, PD2 = x@parameters$PD,
                                                    channelNames = x@parameters$chanNames, incPDofPD = x@parameters$incPDofPD,
                                                    PDofPD = 'pdofpd', removePeptides = x@parameters$removePeps,
                                                    dataType = x@parameters$dataType, modelType = x@parameters$modelType,
                                                    organism = x@parameters$organism,incGeneFile = x@parameters$incGeneFile,
                                                    normaliseData = normalise
              )

             } else {
               x@normData <- normalize_data( dataFrame = x@input,chans = x@parameters$chans,
                                             reps = x@parameters$reps, PD2 = x@parameters$PD,
                                             channelNames = x@parameters$chanNames, incPDofPD = x@parameters$incPDofPD,
                                             PDofPD = 'pdofpd', removePeptides = x@parameters$removePeps,
                                             dataType = x@parameters$dataType, modelType = x@parameters$modelType,
                                             organism = x@parameters$organism,incGeneFile = x@parameters$incGeneFile,
                                             normaliseData = normalise
               )
             }

            return(x)
          }
)
