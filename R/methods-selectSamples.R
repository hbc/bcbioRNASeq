#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase selectSamples
#'
#' @inheritParams general
#'
#' @param transform Apply CPU-intensive DESeq2 transformations. This can
#'   take a long time for large datasets.
#'
#' @return [bcbioRNASeq].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#' load(system.file(
#'     file.path("extdata", "dds.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' selectSamples(bcb, group = "ko")
#'
#' # DESeqDataSet
#' selectSamples(dds, group = "ko")
NULL



# Constructors =================================================================
#' Select Samples Using Dot Arguments
#' @return Character vector of sample identifiers (matching colnames).
#' @noRd
.selectSamples <- function(object, ...) {
    arguments <- list(...)
    invisible(lapply(arguments, assert_is_vector))

    # Match the arguments against the sample metadata
    metadata <- sampleMetadata(object)
    assert_is_data.frame(metadata)

    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        assert_is_subset(column, colnames(metadata))
        argument <- arguments[[a]]
        assert_is_subset(argument, metadata[[column]])
        metadata %>%
            .[.[[column]] %in% argument, "sampleID", drop = TRUE]
    })

    Reduce(f = intersect, x = list)
}



.selectSamples.bcbioRNASeq <- function(  # nolint
    object,
    ...,
    transform = TRUE) {
    samples <- .selectSamples(object, ...)
    object[, samples, transform = transform]
}



.selectSamples.DESeqDataSet <- function(object, ...) {  # nolint
    samples <- .selectSamples(object, ...) %>%
        as.character()
    object <- object[, samples]

    # Relevel the factors in colData
    colData <- colData(object)
    colData[] <- lapply(colData, function(x) {
        if (is.factor(x)) {
            droplevels(x)
        } else {
            x
        }
    })
    colData(object) <- colData

    object
}



# Methods ======================================================================
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioRNASeq"),
    .selectSamples.bcbioRNASeq)



#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("DESeqDataSet"),
    .selectSamples.DESeqDataSet)
