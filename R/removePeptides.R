#' Method to remove peptides from input data of an object of class 'ChemoProtSet'
#'
#' @param x object of class 'ChemoProtSet'
#' @param changePearson option to change the peasrson threshold cut-off parameter
#' @param removePeps boolean value indicating whether peptide removal should take place
#'
#' @return  object of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname removePeptides-methods
#' @examples
#' \dontrun{
#' channelNames <- c("Abundance..F1..126..Control..REP_1",
#'"Abundance..F1..127..Sample..REP_1",  "Abundance..F1..128..Sample..REP_1",
#'  "Abundance..F1..129..Sample..REP_1",  "Abundance..F1..130..Sample..REP_1",
#'"Abundance..F1..131..Sample..REP_1",  "Abundance..F2..126..Control..REP_2",
#' "Abundance..F2..127..Sample..REP_2", "Abundance..F2..128..Sample..REP_2",
#'"Abundance..F2..129..Sample..REP_2",  "Abundance..F2..130..Sample..REP_2",
#'"Abundance..F2..131..Sample..REP_2")
#'ex <- new('ChemoProtSet')
#'ex<- setParameters(x = ex,chansVal = 6, repsVal = 2,
#' dataTypeStr = 'intensity', modelTypeStr = 'linear',
#'  PDBool = FALSE,removePepsBool = FALSE,incPDofPDBool = FALSE,
#'   incGeneFileBool = FALSE,organismStr = 'H.sapiens',
#'    pearsonThrshVal = 0.4)
#'
#'ex<- setData(x = ex, dataFrame = doschedaData,
#'  dataChannels = channelNames,
#'  accessionChannel = "Master.Protein.Accessions",
#'   sequenceChannel = 'Sequence',
#'   qualityChannel = "Qvality.PEP" )
#'ex <- removePeptides(ex,removePeps = FALSE)
#'ex
#'}

setGeneric(name="removePeptides",
           def=function(x, changePearson = NA, removePeps = TRUE)
           {
             standardGeneric("removePeptides")
           }
)
#' @rdname removePeptides-methods
#' @aliases removePeptides,ANY,ANY-method
#'
setMethod(f="removePeptides",
          signature="ChemoProtSet",
          definition=function(x, changePearson = NA, removePeps = TRUE)
          {
            if(!is.na(changePearson)){
              x@parameters$pearsonThrsh <- changePearson
            }

            if(removePeps == FALSE){
              x@parameters$removePeps <- FALSE
            }
            x@normData <- remove_peptides(dataFrame = x@input, chans = x@parameters$chans,
                                                  reps = x@parameters$reps, accessionID = 'Accession',
                                                  chanNames = x@parameters$chanNames, sequenceID = 'Sequence',
                                                  qualityID = 'Quality', incPDofPD = x@parameters$incPDofPD,
                                                  PDofPD = 'pdofpd', removePeptides = x@parameters$removePeps,
                                                  modelType = x@parameters$modelType,
                                                  incGeneFile = x@parameters$incGeneFile, geneFile = x@datasets$GeneID
            )
            return(x)
          }
)
