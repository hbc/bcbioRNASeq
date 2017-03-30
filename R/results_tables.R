#' Differential gene expression results tables
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#' @importFrom basejump setNamesSnake
#'
#' @param bcbio bcbio list object
#' @param results \code{DESeqResults} object
#' @param lfc Log fold change ratio (base 2) cutoff level
#'
#' @return List of DGE results tables
#' @export
#'
#' @examples
#' \dontrun{
#' results_tables(bcbio, res, lfc = 1)
#' }
results_tables <- function(
    bcbio,
    results,
    lfc = 0) {
    if (class(results)[1] != "DESeqResults") {
        stop("results must be a DESeqResults object.")
    }
    name <- deparse(substitute(results))
    alpha <- results@metadata$alpha
    # Running without setting `values = results@rownames` is faster
    annotations <- ensembl_annotations(
        bcbio,
        attributes = c("external_gene_name",
                       "description",
                       "gene_biotype")
    )
    all <- results %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        tibble::as_tibble(.) %>%
        basejump::setNamesSnake(.) %>%
        dplyr::left_join(annotations, by = "ensembl_gene_id") %>%
        dplyr::arrange_(.dots = "padj")
    de <- all %>%
        dplyr::filter_(.dots = ~padj < alpha) %>%
        dplyr::filter_(.dots = ~log2_fold_change < -lfc |
                           log2_fold_change > lfc)
    de_down <- de %>%
        dplyr::filter_(.dots = ~log2_fold_change < 0) %>%
        dplyr::arrange_(.dots = "log2_fold_change")
    de_up <- de %>%
        dplyr::filter_(.dots = ~log2_fold_change > 0) %>%
        # ~-log_fold_change also works.
        # Adding `quote()` here is more readable.
        dplyr::arrange_(.dots = quote(-log2_fold_change))

    writeLines(c(
        name,
        paste("Alpha (FDR) cutoff:", alpha),
        paste("LFC (base 2) cutoff:", lfc),
        paste("Genes evaluated:", nrow(all)),
        # Report non-zero counts?
        paste("Genes upregulated:", nrow(de_up)),
        paste("Genes downregulated:", nrow(de_down))
    ))

    return(list(
        all = all,
        de = de,
        de_down = de_down,
        de_up = de_up,
        dds_name = name
    ))
}
