#' Default boxplot for objects of class ChemoProtSet
#'
#'Description
#'
#' @param x object of class 'ChemoProtSet'
#' @param ... other plotting options
#' @return boxplot for objects of class ChemoProtSet
#' @export
#'
#'@importMethodsFrom  affy boxplot
setMethod(f = "boxplot", signature = "ChemoProtSet", definition = function(x, ...) {
    if (x@parameters$dataType == "intensity") {
        vec <- 1:((x@parameters$chans - 1) * x@parameters$reps)
        
    } else {
        vec <- 1:(x@parameters$chans * x@parameters$reps)
    }
    
    vec <- length(vec)
    palette.bar <- rep(grDevices::terrain.colors(x@parameters$reps), each = ifelse(x@parameters$dataType != 
        "intensity", x@parameters$chans, x@parameters$chans - 1))
    if (nrow(x@normData) == 0) {
        data.merged <- x@input
    } else {
        data.merged <- x@normData
    }
    
    boxplot(data.merged[, 1:vec], col = palette.bar, las = 2, cex.axis = 1, main = c("Box Plots"), 
        ylab = c("Log2(ratios)"), ...)
    legend("topright", legend = paste("rep", 1:x@parameters$reps, sep = " "), fill = grDevices::terrain.colors(x@parameters$reps), 
        cex = 0.6)
    
})
