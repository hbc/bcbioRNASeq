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



# v0.0.25 ======================================================================
#' @rdname deprecated
#' @export
download <- function(...) {
    .Defunct("prepareRNASeqTemplate")
}

#' @rdname deprecated
#' @export
plotGeneDetectionSaturation <- function(...) {
    .Deprecated("plotGeneSaturation")
    plotGeneSaturation(...)
}

#' @rdname deprecated
#' @export
plotDispersion <- function(...) {
    .Deprecated("plotDispEsts")
    plotDispEsts(...)
}



# v0.0.27 ======================================================================
#' @rdname deprecated
#' @export
loadRNASeqRun <- function(...) {
    .Deprecated("loadRNASeq")
    loadRNASeq(...)
}



# v0.1.2 =======================================================================
#' @rdname deprecated
#' @export
plotGeneHeatmap <- function(...) {
    .Deprecated("plotHeatmap")
    plotHeatmap(...)
}



# v0.1.3 =======================================================================
#' @rdname deprecated
#' @export
txi <- function(...) {
    .Defunct()
}



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
#' @importFrom bcbioBase bcbio
#' @export
setMethod(
    "bcbio",
    signature("bcbioRNASeq"),
    function(object, ...) {
        .Defunct()
    }
)

#' @rdname deprecated
#' @importFrom bcbioBase bcbio<-
#' @export
setMethod(
    "bcbio<-",
    signature(
        object = "bcbioRNASeq",
        value = "ANY"
    ),
    function(object, ..., value) {
        .Defunct()
    }
)

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
#' @rdname bcbioRNASeq
#' @usage NULL
#' @export
loadRNASeq <- function(...) {
    .Deprecated("bcbioRNASeq")
    bcbioRNASeq(...)
}



# v0.2.6 =======================================================================
#' @rdname plotMA
#' @usage NULL
#' @export
plotMeanAverage <- function(...) {
    # Soft deprecation, since used in F1000 paper
    plotMA(...)
}



# nocov end
