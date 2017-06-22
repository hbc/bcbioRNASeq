#' YAML utilities
#'
#' @rdname yaml
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param yaml Project summary YAML.
#' @param ... Nested operator keys supplied as dot objects.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @return [DataFrame].
.yaml <- function(yaml, ...) {
    samples <- yaml[["samples"]]
    if (!length(samples)) {
        stop("No sample information in YAML")
    }

    # Check for nested keys, otherwise return NULL
    # TODO Improve the recursion method using sapply in a future update
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
        # Remove legacy duplicate `name` identifier
        nested[["name"]] <- NULL
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
    }
    ) %>%
        rbindlist %>%
        # List can be coerced to data frame using [data.table::rbindlist()] or
        # [dplyr::bind_rows()]. Some YAML files will cause [bind_rows()] to
        # throw `Column XXX can't be converted from integer to character` errors
        # on numeric data, whereas this doesn't happen with [rbindlist()].
        as.data.frame %>%
        remove_na %>%
        arrange(!!sym("description"))
}



#' @rdname yaml
.yaml_metadata <- function(yaml) {
    .yaml(yaml, metadata) %>%
        as.data.frame %>%
        mutate_all(factor) %>%
        mutate(description = as.character(.data[["description"]])) %>%
        as.data.frame %>%
        column_to_rownames("description") %>%
        DataFrame
}



#' @rdname yaml
.yaml_metrics <- function(yaml) {
    metrics <- .yaml(yaml, summary, metrics)
    if (is.null(metrics)) {
        return(NULL)
    }
    characters <- metrics[, c("description",
                              "quality_format",
                              "sequence_length")]
    numerics <- metrics[, setdiff(colnames(metrics), colnames(characters))] %>%
        as.data.frame %>%
        mutate_all(as.numeric)
    bind_cols(characters, numerics) %>%
        as.data.frame %>%
        column_to_rownames("description") %>%
        .[, sort(colnames(.))] %>%
        bind_cols(., characters[, "description", drop = FALSE]) %>%
        DataFrame
}
