#' Sample YAML Metadata Ytilities
#'
#' @rdname yaml
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



#' @rdname yaml
.sample_yaml <- function(yaml, ...) {
    samples <- yaml[["samples"]]
    if (!length(samples)) {
        stop("No sample information in YAML")
    }

    # Check for nested keys, otherwise return NULL
    # Improve recursion method in a future update (lower priority)
    keys <- get_objs_from_dots(dots(...))
    if (!keys[[1L]] %in% names(samples[[1L]])) {
        return(NULL)
    }
    if (length(keys) > 1L) {
        if (!keys[[2L]] %in% names(samples[[1L]][[keys[[1L]]]])) {
            return(NULL)
        }
    }

    lapply(seq_along(samples), function(a) {
        nested <- samples[[a]][[keys]] %>% snake
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
        remove_na %>%
        # Rename `description` to `sample_name`
        rename(sample_name = .data[["description"]]) %>%
        # Set `sample_id` from `sample_name` %>%
        mutate(sample_id = .data[["sample_name"]]) %>%
        .meta_priority_cols
}



#' @rdname yaml
.sample_yaml_metadata <- function(yaml) {
    .sample_yaml(yaml, metadata) %>% .meta_factors
}



#' @rdname yaml
.sample_yaml_metrics <- function(yaml) {
    metrics <- .sample_yaml(yaml, summary, metrics)
    if (is.null(metrics)) {
        return(NULL)
    }
    chr <- metrics %>%
        tidy_select(c(meta_priority_cols,
                      "name",
                      "quality_format",
                      "sequence_length"))
    num <- metrics %>%
        tidy_select(setdiff(colnames(metrics), colnames(chr))) %>%
        mutate_if(is.character, as.numeric)
    bind_cols(chr, num) %>%
        tidy_select(unique(c(meta_priority_cols, sort(colnames(.)))))
}
