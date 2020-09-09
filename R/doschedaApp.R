#' Run shiny application for DOSCHEDA
#'
#' Run a version of the pipline with some extra features and a simple user experience.
#'  The application is documented in detail at here
#' @return Launches shiny application
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom shiny conditionalPanel shinyUI uiOutput sliderInput includeHTML h4 br HTML fluidRow fluidPage radioButtons
#' @import shinydashboard prodlim readxl
#' @export
doschedaApp <- function() {
    appDir <- system.file("shiny-examples", "DOSCHEDA_App", package = "Doscheda")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `Doscheda`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}
