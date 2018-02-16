#' bcbioRNASeq Additional Run Data Accessor (Legacy)
#'
#' Package versions prior to 0.0.27 used `callers` to define the extra bcbio
#' slot. The structure of the object is otherwise the same.
#'
#' @rdname bcbio-legacy
#' @name bcbio-legacy
#' @keywords internal
#'
#' @inherit bcbio
NULL



# Methods ======================================================================
# nocov start

#' @rdname bcbio-legacy
#' @export
setMethod(
    "bcbio",
    signature("bcbioRNADataSet"),
    function(object, type) {
        if (missing(type)) {
            return(slot(object, "callers"))
        }
        if (type %in% names(slot(object, "callers"))) {
            slot(object, "callers")[[type]]
        } else {
            abort(paste(type, "not found"))
        }
    })



#' @rdname bcbio-legacy
#' @export
setMethod(
    "bcbio<-",
    signature(
        object = "bcbioRNADataSet",
        value = "ANY"
    ),
    function(object, type, value) {
        slot(object, "callers")[[type]] <- value
        validObject(object)
        object
    })

# nocov end
