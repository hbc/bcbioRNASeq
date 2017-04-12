#' Get gene annotations from Ensembl
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
#' @param attributes Ensembl attributes. See \code{biomaRt::listAttributes()}.
#' @param filters biomaRt filters. See \code{biomaRt::listFilters()}.
#' @param values Ensembl gene identifier values. Optional but will run faster if
#'   specified.
#'
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_annotations(run, values = "ENSMUSG00000000001")
#' }
ensembl_annotations <- function(
    run,
    attributes = NULL,
    filters = "ensembl_gene_id",
    values = NULL) {
    check_run(run)

    # version = "87"
    ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = paste0(run$organism, "_gene_ensembl"),
        host = "dec2016.archive.ensembl.org"
    )

    # Set biomaRt input defaults
    if (is.null(values)) {
        filters <- ""
        values <- ""
    }

    # Defaults to gene name matching
    if (!is.null(attributes)) {
        attributes <- c("ensembl_gene_id", attributes) %>% unique
    } else {
        attributes <- c("ensembl_gene_id", "external_gene_name")
    }

    annotations <-
        getBM(mart = ensembl,
              attributes = attributes,
              filters = filters,
              values = values) %>%
        arrange_(.dots = colnames(.)) %>%
        set_rownames(.$ensembl_gene_id)

    return(annotations)
}
