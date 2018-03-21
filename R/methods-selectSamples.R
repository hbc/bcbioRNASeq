#' Select Samples
#'
#' @name selectSamples
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase selectSamples
#'
#' @inheritParams general
#' @param transform Apply CPU-intensive DESeq2 transformations. This can
#'   take a long time for large datasets.
#'
#' @return `bcbioRNASeq`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds_small.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq ====
#' x <- selectSamples(bcb_small, day = 7L)
#' show(x)
#' colnames(x)
#'
#' # DESeqDataSet ====
#' x <- selectSamples(dds_small, day = 7L)
#' show(x)
#' colnames(x)
NULL



# Constructors =================================================================
#' Select Samples Using Dot Arguments
#' @return Character vector of sample identifiers (matching colnames).
#' @noRd
.selectSamples <- function(object, ...) {
    arguments <- list(...)
    invisible(lapply(arguments, assert_is_vector))

    # Match the arguments against the sample metadata
    metadata <- colData(object)

    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        assert_is_subset(column, colnames(metadata))
        argument <- arguments[[a]]
        assert_is_subset(argument, metadata[[column]])
        metadata %>%
            .[.[[column]] %in% argument, "sampleID", drop = TRUE]
    })

    Reduce(f = intersect, x = list) %>%
        as.character()
}



# Methods ======================================================================
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioRNASeq"),
    function(object, ..., transform = TRUE) {
        samples <- .selectSamples(object, ...)
        object[, samples, transform = transform]
    }
)



#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("DESeqDataSet"),
    function(object, ...) {
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
)
