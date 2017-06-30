#' Method to fit a model to an object of class 'ChemoProtSet'
#'
#' @param x object of class 'ChemoProtSet'
#'
#' @return  object of class ChemoProtSet
#'
#'
#' @export
#' @docType methods
#' @rdname fitModel-methods
#' @examples
#' channelNames <- c("Abundance..F1..126..Control..REP_1",
#'"Abundance..F1..127..Sample..REP_1",  "Abundance..F1..128..Sample..REP_1",
#'  "Abundance..F1..129..Sample..REP_1",  "Abundance..F1..130..Sample..REP_1",
#'"Abundance..F1..131..Sample..REP_1",  "Abundance..F2..126..Control..REP_2",
#' "Abundance..F2..127..Sample..REP_2", "Abundance..F2..128..Sample..REP_2",
#'"Abundance..F2..129..Sample..REP_2",  "Abundance..F2..130..Sample..REP_2",
#'"Abundance..F2..131..Sample..REP_2")
#' ex <- new('ChemoProtSet')
#' ex<- setParameters(x = ex,chansVal = 6, repsVal = 2,dataTypeStr = 'intensity',
#'  modelTypeStr = 'linear',PDBool = FALSE,removePepsBool = FALSE,
#'  incPDofPDBool = FALSE,incGeneFileBool = FALSE,organismStr = 'H.sapiens', pearsonThrshVal = 0.4)
#' ex<- setData(x = ex, dataFrame = doschedaData, dataChannels = channelNames,
#' accessionChannel = "Master.Protein.Accessions",
#'               sequenceChannel = 'Sequence', qualityChannel = "Qvality.PEP" )
#' ex <- removePeptides(ex,removePeps = FALSE)
#' ex <- runNormalisation(ex)
#' ex <- fitModel(ex)
#' ex
#'
setGeneric(name="fitModel",
           def=function(x)
           {
             standardGeneric("fitModel")
           }
)
#' @rdname fitModel-methods
#'
#'@examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'
#'ex
#' @aliases fitModel,ANY,ANY-method
#'
setMethod(f="fitModel",
          signature="ChemoProtSet",
          definition=function(x)
          {
            if(dim(x@normData)[1] == 0){
              x@finalData <- fit_model(dataFrame = x@input,chans = x@parameters$chans,
                                               reps = x@parameters$reps, PD2 = x@parameters$PD,
                                               sigmoidConc = x@parameters$sigmoidConc, incPDofPD = x@parameters$incPDofPD,
                                               PDofPD = x@parameters$PDofPD,
                                               dataType = x@parameters$dataType, modelType = x@parameters$modelType

              )
            }else{
              x@finalData <- fit_model(dataFrame = x@normData,chans = x@parameters$chans,
                                               reps = x@parameters$reps, PD2 = x@parameters$PD,
                                               sigmoidConc = x@parameters$sigmoidConc, incPDofPD = x@parameters$incPDofPD,
                                               PDofPD = x@parameters$PDofPD,
                                               dataType = x@parameters$dataType, modelType = x@parameters$modelType

              )
            }

            return(x)
          }

)
