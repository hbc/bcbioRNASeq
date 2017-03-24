#' Get gene annotations from Ensembl
#'
#' @author Michael Steinbaugh
#'
#' @keywords internal
#'
#' @import dplyr
#' @importFrom biomaRt getBM listMarts useMart
#'
#' @param organism Organism identifier
#'
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_annotations("mmusculus")
#' }
ensembl_annotations <- function(organism) {
    ensembl_version <- biomaRt::listMarts(host = "useast.ensembl.org") %>%
        dplyr::filter(biomart == "ENSEMBL_MART_ENSEMBL") %>%
        dplyr::select(version) %>% .[[1]]
    message(ensembl_version)

    mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    # mart_options <- biomaRt::listAttributes(mart)

    annotations <-
        biomaRt::getBM(mart = mart,
                       attributes = c("ensembl_gene_id",
                                      "external_gene_name",
                                      "description",
                                      "gene_biotype")) %>%
        dplyr::arrange_(.dots = "ensembl_gene_id") %>%
        set_rownames("ensembl_gene_id")

    save(annotations, file = "data/annotations.rda")

    return(annotations)
}
