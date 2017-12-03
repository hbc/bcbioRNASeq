#' Sample Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @importFrom basejump metrics
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame].
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' metrics(bcb) %>% glimpse()
NULL



# Methods ====
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
