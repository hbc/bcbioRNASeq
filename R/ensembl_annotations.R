#' Get gene annotations from Ensembl
#'
#' @author Michael Steinbaugh
#'
#' @keywords internal
#'
#' @import dplyr
#' @importFrom biomaRt getBM useMart
#'
#' @param bcbio bcbio list object
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
#' ensembl_annotations(bcbio, values = "ENSMUSG00000000001")
#' }
ensembl_annotations <- function(
    bcbio,
    attributes = NULL,
    filters = "ensembl_gene_id",
    values = NULL) {
    check_bcbio_object(bcbio)

    mart <- biomaRt::useMart(
        "ensembl",
        dataset = paste0(bcbio$organism, "_gene_ensembl")
    )
    # attributes <- biomaRt::listAttributes(mart)
    # filters <- biomaRt::listFilters(mart)

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
        biomaRt::getBM(mart = mart,
                       attributes = attributes,
                       filters = filters,
                       values = values) %>%
        dplyr::arrange_(.dots = colnames(.)) %>%
        set_rownames("ensembl_gene_id")

    return(annotations)
}
