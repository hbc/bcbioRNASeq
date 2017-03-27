#' Get gene annotations from Ensembl
#'
#' @author Michael Steinbaugh
#'
#' @keywords internal
#'
#' @import dplyr
#' @importFrom biomaRt getBM useMart
#'
#' @param values Ensembl gene identifier values. Optional but will run faster if
#'   specified.
#' @param filters biomaRt filters. See \code{biomaRt::listFilters()}.
#' @param organism Organism identifier
#'
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_annotations("ENSMUSG00000000001", organism = "mmusculus")
#' }
ensembl_annotations <- function(
    values = NULL,
    filters = "ensembl_gene_id",
    organism) {
    mart <- biomaRt::useMart(
        "ensembl",
        dataset = paste0(organism, "_gene_ensembl")
    )
    # attributes <- biomaRt::listAttributes(mart)
    # filters <- biomaRt::listFilters(mart)

    if (is.null(values)) {
        filters <- ""
        values <- ""
    }

    attributes <- c("ensembl_gene_id",
                    "external_gene_name",
                    "description",
                    "gene_biotype")

    annotations <-
        biomaRt::getBM(mart = mart,
                       attributes = attributes,
                       filters = filters,
                       values = values) %>%
        dplyr::arrange_(.dots = colnames(.)) %>%
        set_rownames("ensembl_gene_id")

    if (dir.exists("data")) {
        save(annotations, file = "data/annotations.rda")
    }

    return(annotations)
}
