#' Show an Object
#'
#' @name show
#' @family S4 Class Definition
#' @author Michael Steinbuagh
#'
#' @inherit methods::show
#'
#' @examples
#' show(bcb_small)
NULL



# Methods ======================================================================
#' @rdname show
#' @export
setMethod(
    "show",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        cat(
            paste(class(object), metadata(object)[["version"]]),
            paste("samples:", ncol(object)),
            paste0(metadata(object)[["level"]], ": ", nrow(object)),
            paste("organism:", metadata(object)[["organism"]]),
            paste("bcbio dir:", metadata(object)[["uploadDir"]]),
            paste("bcbio date:", metadata(object)[["runDate"]]),
            paste("R load date:", metadata(object)[["date"]]),
            sep = "\n"
        )
    }
)
