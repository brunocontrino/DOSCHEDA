% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replicatePlot.R
\docType{methods}
\name{replicatePlot}
\alias{replicatePlot}
\alias{replicatePlot,ChemoProtSet-method}
\alias{replicatePlot,ANY,ANY-method}
\title{Plot replicates between concentrations}
\usage{
replicatePlot(x, conc, repIndex1, repIndex2, ...)

\S4method{replicatePlot}{ChemoProtSet}(x, conc, repIndex1, repIndex2, ...)
}
\arguments{
\item{x}{object of class 'ChemoProtSet'}

\item{conc}{concentration of channel}

\item{repIndex1}{index of replicate on x axis}

\item{repIndex2}{index of replicate on y axis}

\item{...}{options}
}
\value{
Replicate plot for objects of class ChemoProtSet
}
\description{
Plot of Fold Change between replicate i and replicate j at at a given concentration
}
\examples{

ex <- processedExample
ex <- runNormalisation(ex)
ex <- fitModel(ex)
replicatePlot(ex,0,1,2)

}
