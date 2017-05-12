#' Get gene annotations from Ensembl.
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param attributes Ensembl attributes.
#' @param filters biomaRt filters.
#' @param values Ensembl gene identifier values. Optional but will run faster if
#'   specified.
#'
#' @return Data frame.
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_annotations(run, values = "ENSMUSG00000000001")
#' }
#'
#' @seealso
#' \code{\link[biomaRt]{listAttributes}}, \code{\link[biomaRt]{listFilters}}.
ensembl_annotations <- function(
    run,
    attributes = NULL,
    filters = "ensembl_gene_id",
    values = NULL) {
    ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = paste(run$organism, "gene_ensembl", sep = "_"))

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
        arrange(!!!syms(colnames(.))) %>%
        set_rownames(.$ensembl_gene_id)

    return(annotations)
}
