#' Extend SummarizedExperiment Assignment Methods
#'
#' @name SummarizedExperiment
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment assays<-
#'
#' @examples
#' # bcbioRNASeq ====
#' metadata(bcb_small)[["stash"]] <- "test"
#' metadata(bcb_small)[["stash"]]
NULL



# Methods ======================================================================
#' @rdname SummarizedExperiment
#' @export
setMethod(
    "assays<-",
    signature(
        x = "bcbioRNASeq",
        value = "list"
    ),
    getMethod(
        "assays<-",
        signature(
            x = "SummarizedExperiment",
            value = "list"
        )
    )
)



#' @rdname SummarizedExperiment
#' @export
setMethod(
    "assays<-",
    signature(
        x = "bcbioRNASeq",
        value = "SimpleList"
    ),
    getMethod(
        "assays<-",
        signature(
            x = "SummarizedExperiment",
            value = "SimpleList"
        )
    )
)



#' @rdname SummarizedExperiment
#' @export
setMethod(
    "metadata<-",
    signature("bcbioRNASeq"),
    getMethod("metadata<-", "Annotated")
)
