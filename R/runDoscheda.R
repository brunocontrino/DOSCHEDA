#' Wrapper Function to run the entire Doscheda pipeline
#'
#' A wrapper for the whole Doscheda pipeline, if users want to avoid using the separate steps.
#'
#' @param dataFrame data.frame of the input data set
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
#' @param dataChannels column names of dataFrame that correspond to data channels. These should be ordered in the format: rep1_concentration_0, ..., rep1_concentration_n, rep2_concentration_0, ...
#' @param accessionChannel string that is the same as the column name for the protein accessions in dataFrame
#' @param uniquePeps string that is the same as the column name for the number of unique peptides in dataFrame
#' @param sequenceChannel string that is the same as the column name for the peptide sequences in dataFrame
#' @param qualityChannel string that is the same as the column name for the peptide quality score in dataFrame
#' @param pdofpdChannel string that is the same as the column name for the pull-down of pull-down data in dataFrame
#' @param incGeneID boolean value indicating if a protein accession to gene ID file is supplied
#' @param geneIDFile data.frame containing a protein accession to gene ID conversion file
#' @param normType string indicating the type of normalisation that should take place ('loess', 'median', 'none')
#' @return  object of class ChemoProtSet
#'
#' @export
#' @examples
#' channelNames <- c("Abundance..F1..126..Control..REP_1",
#'"Abundance..F1..127..Sample..REP_1",  "Abundance..F1..128..Sample..REP_1",
#'"Abundance..F1..129..Sample..REP_1",  "Abundance..F1..130..Sample..REP_1",
#'"Abundance..F1..131..Sample..REP_1",  "Abundance..F2..126..Control..REP_2",
#'"Abundance..F2..127..Sample..REP_2", "Abundance..F2..128..Sample..REP_2",
#'"Abundance..F2..129..Sample..REP_2",  "Abundance..F2..130..Sample..REP_2",
#' "Abundance..F2..131..Sample..REP_2")
#'
#'ex <- runDoscheda(dataFrame = doschedaData, dataChannels = channelNames,
#' chansVal = 6, repsVal = 2,dataTypeStr = 'intensity',
#' modelTypeStr = 'linear',PDBool = FALSE,removePepsBool = FALSE,
#' accessionChannel = "Master.Protein.Accessions",
#' sequenceChannel = 'Sequence',qualityChannel = "Qvality.PEP",
#' incPDofPDBool = FALSE, incGeneFileBool = FALSE,
#'  organismStr = 'H.sapiens', pearsonThrshVal = 0.4)

   runDoscheda <- function(dataFrame, dataChannels, accessionChannel, chansVal, repsVal, dataTypeStr, modelTypeStr, PDBool = TRUE,
                              removePepsBool = NA, incPDofPDBool = FALSE, PDofPDname = NA,
                              incGeneFileBool = FALSE, organismStr = 'h.sapiens', sigmoidConc = NA,
                              pearsonThrshVal = 0.4,  uniquePeps = NA, sequenceChannel = NA,
                              qualityChannel = NA, pdofpdChannel = NA, incGeneID = FALSE, geneIDFile = NA,
                              normType = 'loess'){

     ex <- new('ChemoProtSet')

     ex<- setParameters(x = ex,chansVal = chansVal, repsVal = repsVal,
                        dataTypeStr = dataTypeStr, modelTypeStr = modelTypeStr,
                        PDBool = PDBool, removePepsBool = removePepsBool,
                        incPDofPDBool = incPDofPDBool,PDofPDname =  PDofPDname,
                        incGeneFileBool = incGeneFileBool, organismStr = organismStr,
                        sigmoidConc = sigmoidConc, pearsonThrshVal = pearsonThrshVal)

     ex<- setData(x = ex, dataFrame = dataFrame, dataChannels = dataChannels,
                  accessionChannel = accessionChannel, uniquePeps = uniquePeps,
                  sequenceChannel =  sequenceChannel, qualityChannel = qualityChannel,
                  pdofpdChannel = pdofpdChannel, incGeneID =  incGeneFileBool,
                  geneIDFile =  geneIDFile)

     if(dataTypeStr  == 'intensity'){
       ex <- removePeptides(ex, removePeps = removePepsBool)
     }

     ex <- runNormalisation(ex, normalise = normType)

     ex <- fitModel(ex)

     return(ex)
   }

