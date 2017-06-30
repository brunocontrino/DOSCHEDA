#' Density plot for objects of class ChemoProtSet
#'
#' Description
#'
#' @param x object of class 'ChemoProtSet'
#' @param rankProteins plot a the set of ranked proteins or plot the density of the channels
#' @param ... other plot options
#'
#' @return density plot for objects of class ChemoProtSet
#'
#' @export
#' @docType methods
#' @rdname densityPlot-methods
#' @examples
#'ex <- processedExample
#'ex <- runNormalisation(ex)
#'ex <- fitModel(ex)
#'densityPlot(ex)

setGeneric(name="densityPlot",
           def=function(x, rankProteins = FALSE, ...)
           {
             standardGeneric("densityPlot")
           }
)

#' @rdname densityPlot-methods
#' @importFrom graphics par legend lines plot points
#' @import grDevices
#' @aliases densityPlot,ANY,ANY-method
#'
setMethod(f= "densityPlot",
          signature="ChemoProtSet",
          definition=function(x, rankProteins = FALSE, ...)
          {
            if(x@parameters$dataType == 'intensity'){
              vec <- 1:((x@parameters$chans -1) * x@parameters$reps)

            }else{
              vec <- 1:(x@parameters$chans * x@parameters$reps)
            }
            pal <- grDevices::rainbow(length(vec))
            leg.nam <- standard_names(x@parameters$chans, x@parameters$reps, x@parameters$dataType)
            if(nrow(x@normData) == 0){
              data.merged <- x@input
            } else {
              data.merged <- x@normData
            }

            missing_val <- 0

            if(rankProteins == TRUE){
              plot(x=rank(data.merged[,1]),y=data.merged[,1], cex.axis=1,
                   main=c(paste("N. of Missing val. ",missing_val ,
                                " \n Change Distribution after LOESS",sep="")),
                   col=pal[1], ylab=c("Log2(ratios)"), xlab = c("Proteins"))
              if(length(vec) > 1){
                for(i in 2:length(vec)){
                  points(x=rank(data.merged[,i]),y=data.merged[,i], cex.axis=1, main=c(""),col=pal[i])
                }
              }
              legend("topleft", legend = leg.nam, fill = pal)

            } else {

              plot(stats::density(x=data.merged[,1]),col= pal[1], main=c(" Density after LOESS"), ylim = c(0,5 ),cex.axis = 1)
              legend("topleft", legend = leg.nam, fill = pal, ...)

              if(length(vec) > 1){
                for(i in 2:length(vec)){
                  lines(stats::density(x=data.merged[,i]), col=pal[i], main = c(""))
                }
              }
            }

          }

)
