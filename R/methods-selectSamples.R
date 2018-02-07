#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase selectSamples
#'
#' @inheritParams AllGenerics
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
#'
#' # bcbioRNASeq
#' selectSamples(bcb, group = "ko")
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
        argument <- arguments[[a]]
        match <- sampleMetadata %>%
            .[.[[column]] %in% argument, "sampleID", drop = TRUE]
        # Check for match failure
        if (!length(match)) {
            warn(paste(
                "Match failure:",
                paste(column, "=", argument)
            ))
            return(NULL)
        }
        match
    })
    samples <- Reduce(f = intersect, x = list)
    if (!length(samples)) {
        warn("No samples matched")
        return(NULL)
    }
    samples %>%
        as.character() %>%
        unique() %>%
        sort()
}



.selectSamples.bcbioRNASeq <- function(  # nolint
    object,
    ...,
    transform = TRUE) {
    samples <- .selectSamples(object, ...)
    object[, samples, transform = transform]
}



#' @importFrom DESeq2 DESeq
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
    DESeq(object)
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
