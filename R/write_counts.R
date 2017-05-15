#' Write gzipped counts.
#'
#' @author Michael Steinbaugh
#'
#' @param ... Count matrices, passed in as dots.
#' @param dir Output directory.
#'
#' @export
write_counts <- function(..., dir = file.path("results", "counts")) {
    names <- as.character(substitute(list(...)))[-1L]
    dots <- list(...)

    # Create the counts output directory
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    message(paste("Writing", toString(names), "to", dir))

    # Iterate across the dots and write CSVs
    lapply(seq_along(names), function(a) {
        name <- names[a]
        counts <- dots[a]
        counts %>%
            as.data.frame %>%
            rownames_to_column("ensembl_gene_id") %>%
            as_tibble %>%
            write_csv(path = file.path(dir, paste0(name, ".csv.gz")))
    }) %>% invisible
}
