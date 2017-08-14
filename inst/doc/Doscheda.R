## ----setParametersEx, eval = TRUE------------------------------------------
library(Doscheda)

 channelNames <- c("Abundance..F1..126..Control..REP_1",
"Abundance..F1..127..Sample..REP_1",  "Abundance..F1..128..Sample..REP_1",
  "Abundance..F1..129..Sample..REP_1",  "Abundance..F1..130..Sample..REP_1",
"Abundance..F1..131..Sample..REP_1",  "Abundance..F2..126..Control..REP_2",
 "Abundance..F2..127..Sample..REP_2", "Abundance..F2..128..Sample..REP_2",
"Abundance..F2..129..Sample..REP_2",  "Abundance..F2..130..Sample..REP_2",
"Abundance..F2..131..Sample..REP_2")
ex <- new('ChemoProtSet')
ex<- setParameters(x = ex,chansVal = 6, repsVal = 2,
dataTypeStr = 'intensity', modelTypeStr = 'linear',
  PDBool = FALSE,removePepsBool = FALSE,incPDofPDBool = FALSE,
   incGeneFileBool = FALSE,organismStr = 'H.sapiens',
    pearsonThrshVal = 0.4)

## ----setDataEx, eval = TRUE------------------------------------------------

ex<- setData(x = ex, dataFrame = Doscheda::doschedaData,
  dataChannels = channelNames,
  accessionChannel = "Master.Protein.Accessions",
   sequenceChannel = 'Sequence',
   qualityChannel = "Qvality.PEP" )


## ----peptideRemovalExample, eval = TRUE------------------------------------

ex <- removePeptides(ex,removePeps = FALSE)


## ----fitModelEx, eval = TRUE, include=FALSE--------------------------------
ex<- runNormalisation(ex)
ex <- fitModel(ex)

## ----runDosEx, eval = FALSE------------------------------------------------
#  
#   channelNames <- c("Abundance..F1..126..Control..REP_1",
#  "Abundance..F1..127..Sample..REP_1",  "Abundance..F1..128..Sample..REP_1",
#  "Abundance..F1..129..Sample..REP_1",  "Abundance..F1..130..Sample..REP_1",
#  "Abundance..F1..131..Sample..REP_1",  "Abundance..F2..126..Control..REP_2",
#  "Abundance..F2..127..Sample..REP_2", "Abundance..F2..128..Sample..REP_2",
#  "Abundance..F2..129..Sample..REP_2",  "Abundance..F2..130..Sample..REP_2",
#   "Abundance..F2..131..Sample..REP_2")
#  
#  ex <- runDoscheda(dataFrame = doschedaData, dataChannels = channelNames,
#   chansVal = 6, repsVal = 2,dataTypeStr = 'intensity',
#   modelTypeStr = 'linear',PDBool = FALSE,removePepsBool = FALSE,
#   accessionChannel = "Master.Protein.Accessions",
#   sequenceChannel = 'Sequence',qualityChannel = "Qvality.PEP",
#   incPDofPDBool = FALSE, incGeneFileBool = FALSE,
#    organismStr = 'H.sapiens', pearsonThrshVal = 0.4)
#  
#  runDoscheda()
#  

## ----plot1-----------------------------------------------------------------
plot(ex)

## ----boxPlot---------------------------------------------------------------
boxplot(ex)

## ----corrplot--------------------------------------------------------------
corrPlot(ex)

## ----pcaplot---------------------------------------------------------------
pcaPlot(ex)

## ----repplot---------------------------------------------------------------
replicatePlot(ex,conc = 0, repIndex1 = 1, repIndex2 = 2)

## ----colcanolot------------------------------------------------------------
volcanoPlot(ex)

## ----meanolcanolot---------------------------------------------------------
meanSdPlot(ex)

## ----exampleR, eval = FALSE------------------------------------------------
#   channelNames <- c("Abundance..F1..126..Control..REP_1",
#  "Abundance..F1..127..Sample..REP_1",  "Abundance..F1..128..Sample..REP_1",
#    "Abundance..F1..129..Sample..REP_1",  "Abundance..F1..130..Sample..REP_1",
#  "Abundance..F1..131..Sample..REP_1",  "Abundance..F2..126..Control..REP_2",
#   "Abundance..F2..127..Sample..REP_2", "Abundance..F2..128..Sample..REP_2",
#  "Abundance..F2..129..Sample..REP_2",  "Abundance..F2..130..Sample..REP_2",
#  "Abundance..F2..131..Sample..REP_2")
#   ex <- new('ChemoProtSet')
#   ex<- setParameters(x = ex,chansVal = 6, repsVal = 2,dataTypeStr = 'intensity',
#    modelTypeStr = 'linear',PDBool = FALSE,removePepsBool = FALSE,
#    incPDofPDBool = FALSE,incGeneFileBool = FALSE,organismStr = 'H.sapiens', pearsonThrshVal = 0.4)
#   ex<- setData(x = ex, dataFrame = doschedaData, dataChannels = channelNames,
#   accessionChannel = "Master.Protein.Accessions",
#                 sequenceChannel = 'Sequence', qualityChannel = "Qvality.PEP" )
#  ex <- removePeptides(ex,removePeps = FALSE)
#  ex <- runNormalisation(ex)
#  ex <- fitModel(ex)
#  ex

