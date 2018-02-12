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
#' selectSamples(dds, group = c("ctrl", "ko"))
NULL



# Constructors =================================================================
#' Select Samples Using Dot Arguments
#' @return Character vector of sample identifiers (matching colnames).
#' @noRd
.selectSamples <- function(object, ...) {
    arguments <- list(...)
    checkArguments <- vapply(
        X = arguments,
        FUN = is.vector,
        FUN.VALUE = logical(1L)
    )
    if (!all(isTRUE(as.logical(checkArguments)))) {
        abort("Arguments must be vectors")
    }

    # Match the arguments against the sample metadata
    sampleMetadata <- sampleMetadata(object)
    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        # Check that column is present
        if (!column %in% colnames(sampleMetadata)) {
            abort(paste(column, "isn't present in metadata colnames"))
        }
        argument <- arguments[[a]]
        # Check that all items in argument are present
        if (!all(argument %in% sampleMetadata[[column]])) {
            missing <- argument[which(!argument %in% sampleMetadata[[column]])]
            abort(paste(
                column,
                "metadata column doesn't contain",
                toString(missing)
            ))
        }
        sampleMetadata %>%
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
    samples <- .selectSamples(object, ...)
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
