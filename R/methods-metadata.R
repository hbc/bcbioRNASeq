#' Metadata
#'
#' @name metadata
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @examples
#' # bcbioRNASeq ====
#' metadata(bcb_small)[["stash"]] <- "test"
#' metadata(bcb_small)[["stash"]]
NULL



# Methods ======================================================================
#' @rdname metadata
#' @importFrom S4Vectors metadata<-
setMethod(
    "metadata<-",
    signature("bcbioRNASeq"),
    getMethod("metadata<-", "Annotated")
)
