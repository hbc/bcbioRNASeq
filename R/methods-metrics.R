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



# Methods ======================================================================
#' @rdname metrics
#' @importFrom dplyr left_join mutate_if
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    function(object) {
        metrics <- metadata(object)[["metrics"]] %>%
            mutate_if(is.character, as.factor)
        metadata <- sampleMetadata(object)
        metrics <- left_join(metrics, metadata, by = metadataPriorityCols)
        rownames(metrics) <- metrics[["sampleID"]]
        metrics
    })
