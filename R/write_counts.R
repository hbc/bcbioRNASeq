#' Write gzipped counts to the results directory
#'
#' @author Michael Steinbaugh
#'
#' @param counts Counts matrix
#'
#' @export
write_counts <- function(counts) {
    dir.create(file.path("results", "counts"),
               recursive = TRUE, showWarnings = FALSE)
    name <- deparse(substitute(counts))
    counts %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        write_csv(
            path = file.path("results", "counts",
                             paste0(name, ".csv.gz")))
}
