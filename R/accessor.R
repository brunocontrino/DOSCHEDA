#'Accessor function for the Input
#'
#'Accessor function for the Input slot of a ChemoProtSet object.
#'
#' @param x object of class ChemoProtSet
#'
#' @return  object of class ChemoProtSet
#' @seealso \code{\link{DoschedaSet}}
#' 
#' @examples
#' ex <- new('ChemoProtSet')
#' getInput(ex)
#'
#' @export
#' @docType methods
#' @rdname getInput-methods
#'
setGeneric(name = "getInput", def = function(x) {
    standardGeneric("getInput")
})

#' @rdname getInput-methods
#' @aliases getInput,ANY,ANY-method
setMethod(f = "getInput", signature = "ChemoProtSet", definition = function(x) {
    return(x@input)
})

#'Accessor function for the normData
#'
#'Accessor function for the normData slot of a ChemoProtSet object.
#'
#' @param x object of class ChemoProtSet
#'
#' @return  object of class ChemoProtSet
#' @seealso \code{\link{DoschedaSet}}
#' @examples
#' ex <- new('ChemoProtSet')
#' getNorm(ex)
#'
#' @export
#' @docType methods
#' @rdname getNorm-methods
#'
setGeneric(name = "getNorm", def = function(x) {
    standardGeneric("getNorm")
})

#' @rdname getNorm-methods
#' @aliases getNorm,ANY,ANY-method
setMethod(f = "getNorm", signature = "ChemoProtSet", definition = function(x) {
    return(x@normData)
})

#'Accessor function for the finalData slot.
#'
#'Accessor function for the finalData slot of a ChemoProtSet object.
#'
#' @param x object of class ChemoProtSet
#'
#' @return  object of class ChemoProtSet
#' @seealso \code{\link{DoschedaSet}}
#' @examples
#' ex <- new('ChemoProtSet')
#' getParameters(ex)
#'
#' @export
#' @docType methods
#' @rdname getFinal-methods
#'
setGeneric(name = "getFinal", def = function(x) {
    standardGeneric("getFinal")
})

#' @rdname getFinal-methods
#' @aliases getFinal,ANY,ANY-method
setMethod(f = "getFinal", signature = "ChemoProtSet", definition = function(x) {
    return(x@finalData)
})

#'Accessor function for the parameters slot.
#'
#'Accessor function for the parameters slot of a ChemoProtSet object.
#'
#' @param x object of class ChemoProtSet
#'
#' @return  object of class ChemoProtSet
#' @seealso \code{\link{DoschedaSet}}
#' @examples
#' ex <- new('ChemoProtSet')
#' getParameters(ex)
#'
#' @export
#' @docType methods
#' @rdname getParameters-methods
#'
setGeneric(name = "getParameters", def = function(x) {
    standardGeneric("getParameters")
})

#' @rdname getParameters-methods
#' @aliases getParameters,ANY,ANY-method
setMethod(f = "getParameters", signature = "ChemoProtSet", definition = function(x) {
    return(x@parameters)
})

#'Accessor function for the datasets slot.
#'
#'Accessor function for the datasets slot of a ChemoProtSet object.
#'
#' @param x object of class ChemoProtSet
#'
#' @return  object of class ChemoProtSet
#' @seealso \code{\link{DoschedaSet}}
#' @examples
#' ex <- new('ChemoProtSet')
#' getDatasets(ex)
#'
#' @export
#' @docType methods
#' @rdname getDatasets-methods
#'
setGeneric(name = "getDatasets", def = function(x) {
    standardGeneric("getDatasets")
})

#' @rdname getDatasets-methods
#' @aliases getDatasets,ANY,ANY-method
setMethod(f = "getDatasets", signature = "ChemoProtSet", definition = function(x) {
    return(x@datasets)
})
