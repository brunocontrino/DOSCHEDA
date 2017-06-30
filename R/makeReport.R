#' Create report from 'ChemProtSet' object
#'
#' Genereate a report that includes several plots and descriptions for an experiment that has been anaylsed using Doscheda
#'
#' @param x Object of class 'ChemoProtSet'
#'
#' @return html report of processed `ChemoProtSet` object
#' @export
#' @examples
#'\dontrun{
#' ex<- new('ChemoProtSet')
#' makeReport(ex)
#'}
makeReport <- function(x){

  temp <- system.file("rmd/report.Rmd", package="Doscheda")
  params <- list(object = x)
  rmarkdown::render(temp, output_file = "report.html", output_dir = '.',
                    params = params,
                    envir = new.env())

}
