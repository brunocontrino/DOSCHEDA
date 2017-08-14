#'Volcano plot for objects of class ChemoProtSet
#'
#'Volcano plots designed to be run on objects of class `ChemoProtSet` when a linear model has been applied.
#'
#' @param x object of class 'ChemoProtSet'
#' @param coefficient coefficient of linear model to be plotted ('slope','intercept','quadratic')
#' @param avExprs average expression cutoff
#' @param pVal p-value cut-off
#' @param ... other plotting options
#'
#' @return volcano plot for objects of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname volcanoPlot-methods
#' @seealso \code{\link{DoschedaSet}}
#' @examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'volcanoPlot(ex)
#'

setGeneric(name = "volcanoPlot", def = function(x, coefficient = "slope", avExprs = 0.2, pVal = 0.05, 
    ...) {
    standardGeneric("volcanoPlot")
})

#' @rdname volcanoPlot-methods
#' @import calibrate
#' @importFrom graphics plot.new
#' @aliases volcanoPlot,ANY,ANY-method
#'
setMethod(f = "volcanoPlot", signature = "ChemoProtSet", definition = function(x, coefficient = "slope", 
    avExprs = 0.2, pVal = 0.05, ...) {
    res <- x@finalData
    # avgthr=0.2 #sign threshold for the averege fold change 0.3(log2) is 1.3 FC
    
    P.Value_intercept <- NULL
    P.Value_quadratic <- NULL
    P.Value_slope <- NULL
    
    if (coefficient == "slope") {
        # Make a basic volcano plot
        with(res, plot(res$AveExpr, -log10(res$P.Value_slope), pch = 20, main = "Volcano plot (slope pval )", 
            xlab = c("Log2_AvgFC"), ylab = c("-Log10(Pval)"), xlim = c(-abs(max(res$AveExpr) + 1), 
                abs(max(res$AveExpr) + 1))))
        
        # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
        s = subset(res, P.Value_slope < pVal)
        with(s, points(s$AveExpr, -log10(s$P.Value_slope), pch = 20, col = "red"))
        
        s = subset(res, abs(res$AveExpr) > avExprs)
        with(s, points(s$AveExpr, -log10(s$P.Value_slope), pch = 20, col = "orange"))
        
        s = subset(res, P.Value_slope < pVal & abs(res$AveExpr) > avExprs)
        with(s, points(s$AveExpr, -log10(s$P.Value_slope), pch = 20, col = "green"))
        # Label points with the textxy function from the calibrate plot
        s = subset(res, P.Value_slope < pVal & abs(res$AveExpr) > avExprs)
        with(s, textxy(s$AveExpr, -log10(s$P.Value_slope), labs = s$GeneID, cex = 1))
    } else if (coefficient == "intercept") {
        
        #### intercept
        
        with(res, plot(res$AveExpr, -log10(res$P.Value_intercept), pch = 20, main = "Volcano plot (Intercept pval )", 
            xlab = c("Log2_AvgFC"), ylab = c("-Log10(Pval)"), xlim = c(-abs(max(res$AveExpr) + 1), 
                abs(max(res$AveExpr) + 1))))
        
        # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
        s = subset(res, P.Value_intercept < pVal)
        with(s, points(s$AveExpr, -log10(s$P.Value_intercept), pch = 20, col = "red"))
        
        s = subset(res, abs(res$AveExpr) > avExprs)
        with(s, points(s$AveExpr, -log10(s$P.Value_intercept), pch = 20, col = "orange"))
        
        s = subset(res, P.Value_intercept < pVal & abs(res$AveExpr) > avExprs)
        with(s, points(s$AveExpr, -log10(s$P.Value_intercept), pch = 20, col = "green"))
        
        # Label points with the textxy function from the calibrate plot
        s = subset(res, P.Value_intercept < pVal & abs(res$AveExpr) > avExprs)
        with(s, textxy(s$AveExpr, -log10(s$P.Value_intercept), labs = s$GeneID, cex = 1))
    } else if (coefficient == "quadratic") {
        
        #### quad
        
        with(res, plot(res$AveExpr, -log10(res$P.Value_quadratic), pch = 20, main = "Volcano plot (Quadratic pval )", 
            xlab = c("Log2_AvgFC"), ylab = c("-Log10(Pval)"), xlim = c(-abs(max(res$AveExpr) + 1), 
                abs(max(res$AveExpr) + 1))))
        
        # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
        s = subset(res, P.Value_quadratic < pVal)
        with(s, points(s$AveExpr, -log10(s$P.Value_quadratic), pch = 20, col = "red"))
        
        s = subset(res, abs(res$AveExpr) > avExprs)
        with(s, points(s$AveExpr, -log10(s$P.Value_quadratic), pch = 20, col = "orange"))
        
        s = subset(res, P.Value_quadratic < pVal & abs(res$AveExpr) > avExprs)
        with(s, points(s$AveExpr, -log10(s$P.Value_quadratic), pch = 20, col = "green"))
        
        # Label points with the textxy function from the calibrate plot
        s = subset(res, P.Value_quadratic < pVal & abs(res$AveExpr) > avExprs)
        with(s, textxy(s$AveExpr, -log10(s$P.Value_quadratic), labs = s$GeneID, cex = 0.9))
    } else {
        message("coefficient choice not accepted. Please choose from: intercept, slope, quadratic")
    }
    legend("bottomleft", title = "Legend", cex = 0.7, c("Not significant", paste("P.Value", pVal, 
        sep = " "), paste("AvgFC >", avExprs, sep = ""), paste("P.Value", pVal, "& AvgFC >", avExprs, 
        sep = "")), col = c("black", "red", "orange", "green"), horiz = FALSE, pch = c(19))
    
})
