#' Extend SummarizedExperiment Assignment Methods
#'
#' @name SummarizedExperiment
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment assays<- colData<-
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
    "colData<-",
    signature(
        x = "bcbioRNASeq",
        value = "DataFrame"
    ),
    getMethod(
        "colData<-",
        signature(
            x = "SummarizedExperiment",
            value = "DataFrame"
        )
    )
)



#' @rdname SummarizedExperiment
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "bcbioRNASeq",
        value = "character"
    ),
    getMethod(
        "interestingGroups<-",
        signature(
            object = "SummarizedExperiment",
            value = "character"
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



#' @rdname SummarizedExperiment
#' @export
setMethod(
    "[[<-",
    signature(
        x = "bcbioRNASeq",
        i = "ANY",
        j = "missing",
        value = "ANY"
    ),
    getMethod(
        "[[<-",
        signature(
            x = "SummarizedExperiment",
            i = "ANY",
            j = "missing",
            value = "ANY"
        )
    )
)
