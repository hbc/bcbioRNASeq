#' Sample YAML Metadata Ytilities
#'
#' @rdname internal-yaml
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param yaml Project summary YAML.
#' @param ... Nested operator keys supplied as dot objects.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @return [tibble].



#' @rdname internal-yaml
.sampleYAML <- function(yaml, ...) {
    samples <- yaml[["samples"]]
    if (!length(samples)) {
        stop("No sample information in YAML")
    }

    # Check for nested keys, otherwise return NULL
    # Improve recursion method in a future update (lower priority)
    keys <- dots(..., character = TRUE)
    if (!keys[[1L]] %in% names(samples[[1L]])) {
        return(NULL)
    }
    if (length(keys) > 1L) {
        if (!keys[[2L]] %in% names(samples[[1L]][[keys[[1L]]]])) {
            return(NULL)
        }
    }

    lapply(seq_along(samples), function(a) {
        nested <- samples[[a]][[keys]]
        # Set the description
        nested[["description"]] <- samples[[a]][["description"]]
        if (rev(keys)[[1L]] == "metadata") {
            if (is.null(nested[["batch"]])) {
                nested[["batch"]] <- NA
            }
            if (length(nested[["phenotype"]])) {
                if (grepl("^$", nested[["phenotype"]])) {
                    nested[["phenotype"]] <- NA
                }
            }
        }
        nested
    }) %>%
        # List can be coerced to data frame using [data.table::rbindlist()] or
        # [dplyr::bind_rows()]. Some YAML files will cause [bind_rows()] to
        # throw `Column XXX can't be converted from integer to character` errors
        # on numeric data, whereas this doesn't happen with [rbindlist()].
        rbindlist %>%
        as("tibble") %>%
        camel %>%
        removeNA %>%
        # Rename `description` to `sampleName`
        rename(sampleName = .data[["description"]]) %>%
        # Set `sampleID` from `sampleName` %>%
        mutate(sampleID = .data[["sampleName"]]) %>%
        .metaPriorityCols
}



#' @rdname internal-yaml
.sampleYAMLMetadata <- function(yaml) {
    .sampleYAML(yaml, metadata) %>% .metaFactors
}



#' @rdname internal-yaml
.sampleYAMLMetrics <- function(yaml) {
    metrics <- .sampleYAML(yaml, summary, metrics)
    if (is.null(metrics)) {
        return(NULL)
    }
    chr <- metrics %>%
        .[, c(metaPriorityCols, "name", "qualityFormat", "sequenceLength")]
    num <- metrics %>%
        .[, setdiff(colnames(metrics), colnames(chr))] %>%
        mutate_if(is.character, as.numeric)
    bind_cols(chr, num) %>%
        .[, unique(c(metaPriorityCols, sort(colnames(.))))]
}
