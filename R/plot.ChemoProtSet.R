#' Default plot for objects of class ChemoProtSet
#'
#'Description
#'
#' @param x object of class 'ChemoProtSet'
#' @param sigmoidCoef the sigmoidal coeffcient, one of ('difference', 'slope', 'rb50'). Obselete if modelType is 'linear'
#' @param ... other plotting options
#' @return plot for objects of class ChemoProtSet
#' @import ggplot2
#' @import gridExtra
#' @import reshape2
#' @export

plot.ChemoProtSet <- function(x, sigmoidCoef = 'rb50', ...) {

  inherits(x, 'ChemoProtSet')

  if(x@parameters$modelType == 'linear'){
    data.merged<- x@finalData
    ..count.. <- NULL
    m0 <- ggplot2::ggplot(data.merged, ggplot2::aes(x=data.merged$P.Value_slope))
    m0<-m0 + geom_histogram(ggplot2::aes(fill = ..count..),binwidth = 0.01) +
      ggplot2::scale_fill_gradient("Count", low = "green", high = "red")+
      ggplot2::xlab("P.val slope")

    ..count.. <- NULL
    m1 <- ggplot2::ggplot(data.merged, ggplot2::aes(x=data.merged$P.Value_intercept))
    m1<- m1 + geom_histogram(ggplot2::aes(fill = ..count..),binwidth = 0.01) +
      ggplot2::scale_fill_gradient("Count", low = "green", high = "red")+
      ggplot2::xlab("Pval intercept")

    m2 <- ggplot2::ggplot(data.merged, ggplot2::aes(x=data.merged$P.Value_quadratic))
    m2 <- m2 + geom_histogram(ggplot2::aes(fill = ..count..),binwidth = 0.01) +
      ggplot2::scale_fill_gradient("Count", low = "green", high = "red")+
      ggplot2::xlab("Pval quadratic")

    gridExtra::grid.arrange(m0,m1,m2)

  } else {



    data_merged_2 <- x@finalData
    conc<- x@parameters$sigmoidConc
    topperc <- 30 #difference in % between top and bottom

    if(x@parameters$dataType == 'intensity'){

      pred.names <- paste0('predX',1:(x@parameters$chans -1))
      final.Names <- paste0('rep1_C',0:(x@parameters$chans - 2))
      diffinter<- data_merged_2[(data_merged_2$predX1 - data_merged_2[,paste("predX",(x@parameters$chans-1),sep = "")]) > topperc & data_merged_2$predX1 <= 100, ]

    } else {
      pred.names <- paste0('predX',1:(x@parameters$chans))
      final.Names <- paste0('rep1_C',0:(x@parameters$chans - 1))
      diffinter<- data_merged_2[(data_merged_2$predX1 - data_merged_2[,paste("predX",(x@parameters$chans),sep = "")]) > topperc & data_merged_2$predX1 <= 100, ]

    }


    if(sigmoidCoef == 'difference'){

      if(nrow(diffinter) > 0){
        Diff_Top_bottom_pred <- shape_for_ggplot_pred(diffinter,log2(conc),pred.names)
        Diff_Top_bottom_perc <- shape_for_ggplot_perc(diffinter,log2(conc),final.Names)
        what<-c("(Top - Bottom) >")
        GeneID <- factor(Diff_Top_bottom_pred$GeneID)
        value<- NULL
        Diff_Top_bottom<-ggplot2::ggplot()+
          ggplot2::geom_line(data = Diff_Top_bottom_pred, ggplot2::aes(x=x,y=value, colour = GeneID), size = 1) +
          ggplot2::geom_point(data = Diff_Top_bottom_perc, ggplot2::aes(x=x,y=value, colour=Diff_Top_bottom_perc$GeneID)) +
          ggplot2::labs(title=paste(what,topperc,sep=""))

        Diff_Top_bottom
      }else{
        Diff_Top_bottom<-ggplot2::ggplot()+
          ggplot2::labs(title=paste("No significant Top-Bottom >" ,topperc,"%","\n","has been found", sep=""))
      }
    } else if(sigmoidCoef == 'slope'){

      ## next plot (SLOPE)

      top<-15 #max prot to plot


      if(x@parameters$dataType == 'intensity'){
        pred.names <- paste0('predX',1:(x@parameters$chans -1))
        final.Names <- paste0('rep1_C',0:(x@parameters$chans - 2))
      }else{
        pred.names <- paste0('predX',1:(x@parameters$chans))
        final.Names <- paste0('rep1_C',0:(x@parameters$chans - 1))
      }

      data_merged_2 <- x@finalData

      #Here make the subselections for using the ggplot functions SLOPE
      slope<-stats::na.omit(data_merged_2[data_merged_2$SlopePval<0.05 ,])
      slope_ordered<-stats::na.omit(slope[order(slope$SlopePval, decreasing = FALSE),][1:top,])
      if(nrow(slope_ordered)>0){
        slope_pred<-shape_for_ggplot_pred(slope_ordered,log10(conc),pred.names)
        slope_perc<- shape_for_ggplot_perc(slope_ordered,log10(conc),final.Names)
        what<-c("Slope (p.val) ")
        GeneID <- factor(slope_pred$GeneID)
        Slope_pl<-ggplot2::ggplot()+
          ggplot2::geom_line(data = slope_pred, ggplot2::aes(x=x,y = value, colour = GeneID), size = 1) +
          ggplot2::geom_point(data = slope_perc, ggplot2::aes(x=x,y = value,colour=slope_perc$GeneID))+
          ggplot2::labs(title=paste(what,"Top",top,sep=""))
        Slope_pl
      }else{Slope_pl<-ggplot2::ggplot()+
        ggplot2::labs(title="No significant Sigmoidal Slope has been found")
      }

    } else if(sigmoidCoef == 'rb50'){

      ##### next plot (RB50)


      top<-15 #max prot to plot

      # RB50<-na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval<0.05 & data_merged_2$predX1-data_merged_2$predX9 >0 & data_merged_2$predX1 <= 100,])
      RB50 <- data.frame(stats::na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval < 0.05
                                                      & data_merged_2$predX1-data_merged_2[,paste0('predX',(x@parameters$chans - 1))] >0 & data_merged_2$predX1 <= 100,]))


      RB50_ordered<- stats::na.omit(RB50[order(RB50$RB50Pval, decreasing = FALSE),][1:top,])

      if(nrow(RB50_ordered)>0){
        RB50_pred<-shape_for_ggplot_pred(RB50_ordered,log10(conc),pred.names)
        RB50_perc<-shape_for_ggplot_perc(RB50_ordered,log10(conc),final.Names)
        what<-c("RB50 (p.val) ")
        GeneID <- factor(RB50_pred$GeneID)
        RB50_pl<-ggplot2::ggplot()+
          ggplot2::geom_line(data = RB50_pred, ggplot2::aes(x=x,y=value, colour = GeneID), size = 1) +
          ggplot2::geom_point(data = RB50_perc, ggplot2::aes(x=x,y=value,colour=RB50_perc$GeneID))+
          ggplot2::labs(title=paste(what,"Top",top,sep=""))
        RB50_pl
      }else{
        RB50_pl<-ggplot2::ggplot()+
          ggplot2::labs(title="No significant RB50 has been found")
        print(RB50_pl)
      }
    }else{
      message('sigmoidCoef not accepted please enter one of: "difference", "slope" or "rb50"')
    }



  }

}
