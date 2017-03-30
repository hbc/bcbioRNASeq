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
#' @param results \code{DESeqResults} or \code{DESeqDataSet} object
#' @param alpha Alpha cutoff level
#' @param lfc Log fold change ratio (base 2) cutoff level
#' @param organism Organism
#' @param print Print summary
#'
#' @return List of DGE results tables
#' @export
#'
#' @examples
#' \dontrun{
#' results_tables(res, alpha = 0.05, lfc = 1, organism = "hsapiens")
#' }
results_tables <- function(
    bcbio,
    results,
    alpha = 0.05,
    lfc = 0,
    print = TRUE) {
    # Add test for alpha as numeric, less than 1
    # Add test for lfc as numeric, non-negative

    # Coerce `DESeqDataSet` to `DESeqResults` if necessary
    if (class(results)[1] == "DESeqDataSet") {
        results <- DESeq2::results(results, alpha = alpha)
    }
    if (class(results)[1] != "DESeqResults") {
        stop("results must be a DESeqResults object.")
    }

    annotations <- ensembl_annotations(
        bcbio,
        attributes = c("external_gene_name",
                       "description",
                       "gene_biotype"),
        values = results@rownames
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

    if (isTRUE(print)) {
        writeLines(c(
            paste("Alpha (FDR) cutoff:", alpha),
            paste("LFC (base 2) cutoff:", lfc),
            # Add non-zero count here?
            paste("Genes evaluated:", nrow(all)),
            paste("Genes upregulated:", nrow(de_up)),
            paste("Genes downregulated:", nrow(de_down))
        ))
    }

    return(list(
        all = all,
        de = de,
        de_down = de_down,
        de_up = de_up
    ))
}
