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
#' @importFrom dplyr left_join mutate_if
#' @importFrom S4Vectors metadata
.metrics.bcbioRNASeq <- function(object) {  # nolint
    data <- metadata(object)[["metrics"]]
    assert_is_data.frame(data)
    data <- mutate_if(data, is.character, as.factor)
    metadata <- sampleMetadata(object)
    assert_is_data.frame(metadata)
    data <- left_join(data, metadata, by = metadataPriorityCols)
    rownames(data) <- data[["sampleID"]]
    data
}



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    .metrics.bcbioRNASeq)
