#' Method to set parameters for a ChemoProtSet
#'
#' Give the ChemoProtSet object the correct parameters for a given experiment in order to successfully run the pipelin
#'
#' @param x object of class 'ChemoProtSet'
#' @param chansVal number of channels / concentrations in experiment
#' @param repsVal number of replicates in experiment
#' @param dataTypeStr string describing the data type of input data set. This can be 'LFC' for log fold-changes, 'FC' for fold-changes and 'intensity' for peptide intensities
#' @param modelTypeStr string describing the type of model applied. This can be 'linear' for a linear model or 'sigmoid' for a sigmoidal model
#' @param PDBool boolean value indicating if the input data is from Proteome Discoverer 2.1 or not
#' @param removePepsBool boolean value indicating if peptide removal will take place. Only valid if input data is peptide intensities
#' @param incPDofPDBool boolean value indicating if the input data contais a pull-down of pull-down colum
#' @param PDofPDname string with the same name as colulmn containing pull-down of pull-down data. NA if this is not applicable
#' @param incGeneFileBool boolean value indicating if the data requires a protein accession to gene ID conversion file
#' @param organismStr string giving the name of organism. the options are: 'H.sapiens', 'D. melanogaster', 'C. elegans', 'R. norvegicus', 'M. musculus'. This is only needed if PDbool is FALSE
#' @param sigmoidConc vector of numerical values for concentrations of channels in the case of a sigmoidal fit
#' @param pearsonThrshVal numerial value between -1 and 1 which determines the cut-off used to discard peptides during peptide removal
#'
#' @return  object of class ChemoProtSet 
#' @seealso \code{\link{DoschedaSet}}
#'
#' @examples
#' channelNames <- c('Abundance..F1..126..Control..REP_1',
#''Abundance..F1..127..Sample..REP_1',  'Abundance..F1..128..Sample..REP_1',
#''Abundance..F1..129..Sample..REP_1',  'Abundance..F1..130..Sample..REP_1',
#''Abundance..F1..131..Sample..REP_1',  'Abundance..F2..126..Control..REP_2',
#''Abundance..F2..127..Sample..REP_2', 'Abundance..F2..128..Sample..REP_2',
#''Abundance..F2..129..Sample..REP_2',  'Abundance..F2..130..Sample..REP_2',
#' 'Abundance..F2..131..Sample..REP_2')
#'
#'ex <- new('ChemoProtSet')
#'ex<- setParameters(x = ex,chansVal = 6, repsVal = 2,dataTypeStr = 'intensity',
#'modelTypeStr = 'linear',PDBool = FALSE, removePepsBool = FALSE,
#'incPDofPDBool = FALSE, incGeneFileBool = FALSE,
#'organismStr = 'H.sapiens', pearsonThrshVal = 0.4)
#'
#'ex
#' @export
#' @docType methods
#' @rdname setParameters-methods
setGeneric(name = "setParameters", def = function(x, chansVal, repsVal, dataTypeStr, modelTypeStr, 
    PDBool = TRUE, removePepsBool = NA, incPDofPDBool = FALSE, PDofPDname = NA, incGeneFileBool = FALSE, 
    organismStr = "h.sapiens", sigmoidConc = NA, pearsonThrshVal = 0.4) {
    standardGeneric("setParameters")
})

#' @rdname setParameters-methods
#' @aliases setParameters,ANY,ANY-method
#'
setMethod(f = "setParameters", signature = "ChemoProtSet", definition = function(x, chansVal, repsVal, 
    dataTypeStr, modelTypeStr, PDBool = TRUE, removePepsBool = NA, incPDofPDBool = FALSE, PDofPDname = NA, 
    incGeneFileBool = FALSE, organismStr = "h.sapiens", sigmoidConc = NA, pearsonThrshVal = 0.4) {
    
    if (dataTypeStr == "intensity") {
        
        tempNames <- paste("rep", rep(1:repsVal, each = chansVal), "_", rep(paste("C", c("ontrol", 
            0:(chansVal - 2)), sep = ""), repsVal), sep = "")
        
        
        
    } else {
        tempNames <- paste("rep", rep(1:repsVal, each = chansVal), "_", rep(paste("C", 0:(chansVal - 
            1), sep = ""), repsVal), sep = "")
        
    }
    
    x@parameters <- list(chans = chansVal, reps = repsVal, dataType = dataTypeStr, modelType = modelTypeStr, 
        chanNames = tempNames, PD = PDBool, removePeps = removePepsBool, incPDofPD = incPDofPDBool, 
        PDofPD = PDofPDname, incGeneFile = incGeneFileBool, organism = organismStr, sigmoidConc = sigmoidConc, 
        pearsonThrsh = pearsonThrshVal)
    
    return(x)
})
