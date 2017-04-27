#' Read sample 
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#' @param nested_operators Nested operators work as an ordered character vector,
#'     recursing a level down for each entry
#'
#' @export
parse_samples_yaml <- function(run, nested_operators) {
    check_run(run)
    yaml <- run$yaml
    if (is.null(yaml)) {
        stop("Run YAML summary is required")
    }
    samples <- yaml$samples
    if (!length(samples)) {
        stop("No sample information in YAML")
    }
    df <- lapply(seq_along(samples), function(a) {
        # Unlist to named character vector
        c(description = samples[[a]]$description,
          unlist(samples[[a]][[nested_operators]])) %>%
            # Replace empty with NA
            gsub("^$", NA, .) %>%
            # Sanitize names in snake_case
            set_names_snake
    }) %>%
        do.call(rbind, .) %>%
        as.data.frame %>%
        # Order by description
        .[order(.$description), ] %>%
        # Remove columns with only NAs
        .[colSums(!is.na(.)) > 0] %>%
        set_rownames(.$description)
    return(df)
}
