#'Plot showing correlation between all channels across replicates
#'
#' Plot of the correlation between all the channels in the data.
#'
#' @param x object of class 'ChemoProtSet'
#' @param ... corrplot options
#' @return correlation plot for objects of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname corrPlot-methods
#' @examples
#'
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'corrPlot(ex)
#'
#'

setGeneric(name="corrPlot",
           def=function(x, ...)
           {
             standardGeneric("corrPlot")
           }
)

#' @rdname corrPlot-methods
#' @importFrom corrgram corrgram
#' @importFrom stats cor
#' @aliases corrPlot,ANY,ANY-method
#'
setMethod(f= "corrPlot",
          signature="ChemoProtSet",
          definition=function(x, ...)
          {

            if(x@parameters$modelType == 'sigmoid'){
              su <- x@normData

            }else{

              if(x@parameters$reps == 1 ){


                su <- x@normData

              }else{
                su <- x@finalData
              }
            }

            if(x@parameters$dataType == 'intensity'){

              index<- 1:((x@parameters$chans - 1) * x@parameters$reps)
            }else{
              index<- 1:(x@parameters$chans * x@parameters$reps)
            }
            # panel.shadeNtext <- NULL; panel.pie <- NULL;  panel.txt <- NULL ;

            index<- length(index)
            # panel.shadeNtext
            cmat <-cor(su[,1:index],use="pairwise", method = "pearson")
            corrgram(cmat, order=TRUE, lower.panel = function (x, y, corr = NULL, col.regions, ...)
            {
              if (is.null(corr))
                corr <- stats::cor(x, y, use = "
                                   pair")
              ncol <- 14
              pal <- col.regions(ncol)
              col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1,
                                                           length = ncol + 1), include.lowest = TRUE))
              usr <- graphics::par("usr")
              graphics::rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind],
                             border = NA)
              graphics::box(col = "lightgray")
              on.exit(graphics::par(usr))
              graphics::par(usr = c(0, 1, 0, 1))
              r <- formatC(corr, digits = 2, format = "f")
              cex.cor <- .8/graphics::strwidth("-X.xx")
              graphics::text(0.5, 0.5, r, cex = cex.cor)},
                     upper.panel= corrgram::panel.pie, text.panel = corrgram::panel.txt,
                     main="Corrgram Plots" )





          }

)
