#' Sample Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#'
#' @return [data.frame].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' metrics(bcb) %>% glimpse()
NULL



# Constructors =================================================================
.metrics.bcbioRNASeq <- function(object) {  # nolint
    metrics <- metadata(object)[["metrics"]]
    assert_is_data.frame(metrics)
    metrics <- mutate_if(metrics, is.character, as.factor)
    metadata <- sampleMetadata(object)
    metrics <- left_join(metrics, metadata, by = metadataPriorityCols)
    rownames(metrics) <- metrics[["sampleID"]]
    metrics
}



# Methods ======================================================================
#' @rdname metrics
#' @importFrom dplyr left_join mutate_if
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    .metrics.bcbioRNASeq)
