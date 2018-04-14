#' Extend SummarizedExperiment Assignment Methods
#'
#' @name SummarizedExperiment
#' @author Michael Steinbaugh
#' @keywords internal
#' @importFrom SummarizedExperiment assays<- colData<-
#'
#' @examples
#' # bcbioRNASeq ====
#' metadata(bcb_small)[["stash"]] <- "test"
#' metadata(bcb_small)[["stash"]]
NULL



# Methods ======================================================================
#' @rdname SummarizedExperiment
#' @importFrom SummarizedExperiment assays<-
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
#' @importFrom SummarizedExperiment assays<-
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
#' @importFrom SummarizedExperiment colData<-
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
#' @importFrom bcbioBase interestingGroups<-
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
#' @importFrom S4Vectors metadata<-
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
