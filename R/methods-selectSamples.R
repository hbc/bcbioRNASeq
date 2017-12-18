#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#' @author Michael Steinbaugh
#'
#' @importFrom basejump selectSamples
#'
#' @inheritParams AllGenerics
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
    checkCharacter <- vapply(arguments, is.character, FUN.VALUE = logical(1))
    if (!all(isTRUE(as.logical(checkCharacter)))) {
        stop("'Arguments must be character")
    }

    # Match the arguments against the sample metadata
    sampleMetadata <- sampleMetadata(object)
    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        argument <- arguments[[a]]
        match <- sampleMetadata %>%
            .[.[[column]] %in% argument,
              "sampleID",
              drop = TRUE]
        # Check for match failure
        if (!length(match)) {
            warning(paste(
                "Match failure:",
                paste(column, "=", argument)
            ), call. = FALSE)
            return(NULL)
        }
        match
    })
    samples <- Reduce(f = intersect, x = list)
    if (!length(samples)) {
        warning("No samples matched", call. = FALSE)
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
