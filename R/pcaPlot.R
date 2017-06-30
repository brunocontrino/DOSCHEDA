#' PCA of the main data sets contained in a object of class ChemoProtSet
#'
#'Plot of Principal Component Analysis for the first two principal components of the experimental data.
#'
#' @param x object of class 'ChemoProtSet'
#' @param ... other plot options
#'
#' @return PCA plot for objects of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname pcaPlot-methods
#' @examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'pcaPlot(ex)


setGeneric(name="pcaPlot",
           def=function(x, ...)
           {
             standardGeneric("pcaPlot")
           }
)

#' @rdname pcaPlot-methods
#' @import ggplot2
#' @importFrom stats prcomp
#' @aliases pcaPlot,ANY,ANY-method
#' @examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'pcaPlot(ex)
#'
setMethod(f= "pcaPlot",
          signature="ChemoProtSet",
          definition=function(x, ...)
          {
            if( x@parameters$dataType == 'intensity'){

              nchan<- x@parameters$chans - 1

            }else{

              nchan<- x@parameters$chans
            }

            su <- x@finalData

            reps <- x@parameters$reps
            index <- 1:(nchan*reps)
            pca <- prcomp(su[,1:length(index)], scale=FALSE)

            DTA<-data.frame( as.numeric(t(su[,1:length(index)])%*%pca$x[,1]),
                             as.numeric(t(su[,1:length(index)])%*%pca$x[,2]))

            p<-ggplot2::ggplot(DTA, ggplot2::aes(x=DTA$as.numeric.t.su...1.length.index........pca.x...1..,
                                                 y=DTA$as.numeric.t.su...1.length.index........pca.x...2..))

            paste0("PC1", " (", round(pca$sdev[1]/sum(pca$sdev)*100,0), "%)")

            shapeval <- c(15:18,7:12)
            p <- p + ggplot2::geom_point(ggplot2::aes(colour = factor(rep(1:reps,each = (nchan)),labels = paste("Rep",1:reps))[index],
                                                      shape = factor(rep(1:nchan,reps),labels = c(paste("C",0:(nchan-1),sep = "")))[index] ), size = 5 ) +
              scale_shape_manual(values=shapeval[1:nchan]) + ggplot2::labs(x = paste0("PC1", " (", round(pca$sdev[1]/sum(pca$sdev)*100,0), "%)"),
                                                                           y = paste0("PC2", " (", round(pca$sdev[2]/sum(pca$sdev)*100,0), "%)"), title="PCA") + ggplot2::labs(color = "Replicates", shape="Concentration")

            p
          }

)
