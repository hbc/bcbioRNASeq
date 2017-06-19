#' Generate a SummarizedExperiment object
#'
#' @rdname SummarizedExperiment
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param txi tximport list.
#' @param metadata Metadata [SimpleList].
#' @param ensembl Ensembl annotations.
.SummarizedExperiment <- function(txi, metadata, ensembl) {
    colData <- metadata
    rowData <- ensembl[["gene"]] %>%
        .[rownames(txi[["counts"]]), ]
    SummarizedExperiment(
        assays = SimpleList(
            abundance = txi[["abundance"]],
            counts = txi[["counts"]],
            length = txi[["length"]]),
        colData = colData,
        rowData = rowData,
        metadata = SimpleList(
            counts_from_abundance = txi[["countsFromAbundance"]]))
}
