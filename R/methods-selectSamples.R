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
.selectSamples <- function(
    object,
    ...,
    transform = TRUE) {
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
    samples <- sort(unique(samples))
    object[, samples, transform = transform]
}



# Methods ======================================================================
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioRNASeq"),
    .selectSamples)
