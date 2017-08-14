#' Method for attaching and standardising data for objects of class 'ChemoProtSet'
#'
#' This method will subset the orginal data set into the required columns, standardising column names in the process.
#'
#' @param x object of class 'ChemoProtSet'
#' @param dataFrame data.frame of the input data set
#' @param dataChannels column names of dataFrame that correspond to data channels. These should be ordered in the format: rep1_concentration_0, ..., rep1_concentration_n, rep2_concentration_0, ...
#' @param accessionChannel string that is the same as the column name for the protein accessions in dataFrame
#' @param uniquePeps string that is the same as the column name for the number of unique peptides in dataFrame
#' @param sequenceChannel string that is the same as the column name for the peptide sequences in dataFrame
#' @param qualityChannel string that is the same as the column name for the peptide quality score in dataFrame
#' @param pdofpdChannel string that is the same as the column name for the pull-down of pull-down data in dataFrame
#' @param incGeneID boolean value indicating if a protein accession to gene ID file is supplied
#' @param geneIDFile data.frame containing a protein accession to gene ID conversion file
#'
#' @return  object of class ChemoProtSet
#' @seealso \code{\link{DoschedaSet}}
#' @examples
#'
#' channelNames <- c('Abundance..F1..126..Control..REP_1',
#''Abundance..F1..127..Sample..REP_1',  'Abundance..F1..128..Sample..REP_1',
#'  'Abundance..F1..129..Sample..REP_1',  'Abundance..F1..130..Sample..REP_1',
#''Abundance..F1..131..Sample..REP_1',  'Abundance..F2..126..Control..REP_2',
#' 'Abundance..F2..127..Sample..REP_2', 'Abundance..F2..128..Sample..REP_2',
#''Abundance..F2..129..Sample..REP_2',  'Abundance..F2..130..Sample..REP_2',
#' 'Abundance..F2..131..Sample..REP_2')
#'
#'ex <- new('ChemoProtSet')
#'ex<- setParameters(x = ex,chansVal = 6, repsVal = 2,dataTypeStr = 'intensity',
#' modelTypeStr = 'linear',PDBool = FALSE,removePepsBool = FALSE,
#' incPDofPDBool = FALSE,incGeneFileBool = FALSE,organismStr = 'H.sapiens', pearsonThrshVal = 0.4)
#'ex<- setData(x = ex, dataFrame = doschedaData, dataChannels = channelNames,
#' accessionChannel = 'Master.Protein.Accessions',
#' sequenceChannel = 'Sequence',qualityChannel = 'Qvality.PEP')
#'
#'ex
#'
#' @export
#' @docType methods
#' @rdname setData-methods
#'
setGeneric(name = "setData", def = function(x, dataFrame, dataChannels, accessionChannel, uniquePeps = NA, 
    sequenceChannel = NA, qualityChannel = NA, pdofpdChannel = NA, incGeneID = FALSE, geneIDFile = NA) {
    standardGeneric("setData")
})

#' @rdname setData-methods
#' @aliases setData,ANY,ANY-method
setMethod(f = "setData", signature = "ChemoProtSet", definition = function(x, dataFrame, dataChannels, 
    accessionChannel, uniquePeps = NA, sequenceChannel = NA, qualityChannel = NA, pdofpdChannel = NA, 
    incGeneID = FALSE, geneIDFile = NA) {
    
    
    if (x@parameters$dataType == "intensity") {
        if (x@parameters$incPDofPD == TRUE) {
            
            x@input <- dataFrame[, c(accessionChannel, sequenceChannel, qualityChannel, pdofpdChannel, 
                dataChannels)]
            colnames(x@input) <- c("Accession", "Sequence", "Quality", "pdofpd", x@parameters$chanNames)
        } else {
            x@input <- dataFrame[, c(accessionChannel, sequenceChannel, qualityChannel, dataChannels)]
            colnames(x@input) <- c("Accession", "Sequence", "Quality", x@parameters$chanNames)
            
        }
        
    } else {
        
        if (x@parameters$PD == TRUE) {
            uniPep <- colnames(dataFrame)[grep("unique", colnames(dataFrame), ignore.case = TRUE)]
            x@input <- dataFrame[, c(accessionChannel, "Description", uniPep, dataChannels)]
            colnames(x@input) <- c("Accession", "Description", "UniquePeptides", x@parameters$chanNames)
            
        } else {
            
            x@input <- dataFrame[, c(accessionChannel, uniquePeps, dataChannels)]
            colnames(x@input) <- c("Accession", "UniquePeptides", x@parameters$chanNames)
        }
    }
    
    if (incGeneID == TRUE) {
        
        if (nrow(geneIDFile[, grep("Accession", geneIDFile)]) == 0) {
            stop("Your Accession to Gene ID conversion file does not contain a column called 'Accession'. Please ensure it contains a column called 'Accession' and a column called 'GeneID'.")
        } else if (nrow(geneIDFile[, grep("GeneID", geneIDFile)]) == 0) {
            
            stop("Your Accession to Gene ID conversion file does not contain a column called 'GeneID'. Please ensure it contains a column called 'Accession' and a column called 'GeneID'.")
            
        } else {
            
            x@datasets$geneID <- data.frame(Accession = geneIDFile[, "Accession"], GeneID = geneIDFile[, 
                "GeneID"])
            
        }
        
    }
    
    return(x)
})
