#' Plot replicates between concentrations
#'
#' Plot of Fold Change between replicate i and replicate j at at a given concentration
#'
#' @param x object of class 'ChemoProtSet'
#' @param conc concentration of channel
#' @param repIndex1 index of replicate on x axis
#' @param repIndex2 index of replicate on y axis
#' @param ... options
#' @return Replicate plot for objects of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname replicatePlot-methods
#' @examples
#'
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'replicatePlot(ex,0,1,2)
#'
setGeneric(name = "replicatePlot", def = function(x, conc, repIndex1, repIndex2, ...) {
    standardGeneric("replicatePlot")
})

#' @rdname replicatePlot-methods
#' @importFrom graphics text title
#' @aliases replicatePlot,ANY,ANY-method
setMethod(f = "replicatePlot", signature = "ChemoProtSet", definition = function(x, conc, repIndex1, 
    repIndex2, ...) {
    if (x@parameters$reps > 1) {
        
        data.merged <- x@finalData
        
        index <- index_matrix(x@parameters$chans, x@parameters$reps, x@parameters$dataType)
        index <- index[index$concentration == conc, ]
        index <- index[, index$rep1 == repIndex1 & index$rep2 == repIndex2]
        
        val = max(c(max(data.merged[, index$index1], na.rm = TRUE), max(data.merged[, index$index2], 
            na.rm = TRUE)))
        
        
        plot(x = data.merged[, index$index1], y = data.merged[, index$index2], xlim = c(0, val + 
            0.2), col = "green", ylim = c(0, val + 0.2), cex.axis = 1.2, main = c(paste0("LogFC-LogFC Plot: Concentration ", 
            conc)), xlab = paste("rep", index$rep1), ylab = paste("rep", index$rep2))
        text(data.merged[, index$index1], data.merged[, index$index2], labels = data.merged$GeneID, 
            cex = 1, pos = 4)
        lines(x = c(0, val), y = c(0, val), col = "red")
    } else {
        plot.new()
        title(main = "Plot not available: only 1 replicate")
        
    }
    
})
