# nocov start

#' Defunct or Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return No value.
NULL



# v0.2.0 =======================================================================
# annotable deprecation for SummarizedExperiment added to bcbioBase v0.2.0

#' @rdname deprecated
#' @export
meltLog10 <- function(object, ...) {
    .Defunct()
}

#' @rdname deprecated
#' @export
plot53Bias <- function(...) {
    .Deprecated("plot5Prime3PrimeBias")
    plot5Prime3PrimeBias(...)
}

#' @rdname deprecated
#' @export
setMethod(
    "design",
    signature("bcbioRNASeq"),
    function(object, ...) {
        .Defunct()
    }
)

#' @rdname deprecated
#' @importFrom BiocGenerics design<-
#' @export
setMethod(
    "design<-",
    signature(
        object = "bcbioRNASeq",
        value = "formula"
    ),
    function(object, ..., value) {
        .Defunct()
    }
)



# v0.2.2 =======================================================================
#' @rdname deprecated
#' @export
loadRNASeq <- function(...) {
    .Deprecated("bcbioRNASeq")
    bcbioRNASeq(...)
}



# v0.2.6 =======================================================================
#' @rdname deprecated
#' @export
plotMeanAverage <- function(...) {
    # Soft deprecation, since used in F1000 paper.
    plotMA(...)
}



# nocov end
