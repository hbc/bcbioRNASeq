# https://github.com/hbc/CHBUtils/blob/master/R/volcanoPlot.R
# `volcano_density_plot()` requires a `data.frame` with two columns:
# `logFC` and `Adjusted.Pvalue`


#' Volcano plot
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#' @importFrom CHBUtils volcano_density_plot
#'
#' @param bcbio bcbio list object
#' @param dds DESeq2 data set
#' @param alpha Alpha level cutoff for coloring
#' @param lfc Log fold change ratio (base 2) cutoff for coloring
#' @param text_alpha Alpha level cutoff for text labels
#' @param text_lfc Log fold change ratio (base 2) cutoff for text labels
#'
#' @return Volcano plot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_volcano(bcbio,
#'              dds,
#'              alpha = 0.05,
#'              lfc = 1)
#' }
plot_volcano <- function(bcbio,
                         dds,
                         alpha = 0.05,
                         lfc = 1,
                         text_alpha = 1e-8,
                         text_lfc = 1) {
    check_bcbio_object(bcbio)
    # We can modify this function to support other methods in the future.
    # For now, it supports DESeq2 results.
    if (class(dds)[1] != "DESeqDataSet") {
        stop("A DESeqDataSet object is required.")
    }

    # Prepare data frame for `CHBUtils::volcano_density_plot()`. Note that
    # DESeq2 result table columns must be renamed.
    df <- dds %>%
        DESeq2::results(.) %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        # Filter zero counts, select columns
        dplyr::filter_(.dots = ~!is.na(log2FoldChange)) %>%
        dplyr::select_(.dots = c("ensembl_gene_id",
                                 "log2FoldChange",
                                 "padj")) %>%
        dplyr::rename_(.dots = c(
            "Adjusted.Pvalue" = "padj",
            "logFC" = "log2FoldChange"
        )) %>%
        set_rownames("ensembl_gene_id") %>%
        dplyr::select_(.dots = c("logFC", "Adjusted.Pvalue"))

    plot_text <- df %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        dplyr::filter_(.dots = ~Adjusted.Pvalue < text_alpha) %>%
        dplyr::filter_(.dots = ~logFC < -text_lfc | logFC > text_lfc) %>%
        dplyr::left_join(
            ensembl_annotations(bcbio, values = .$ensembl_gene_id),
            by = "ensembl_gene_id"
        ) %>%
        dplyr::rename_(.dots = c("name" = "external_gene_name")) %>%
        dplyr::select_(.dots = c("ensembl_gene_id",
                                 "logFC",
                                 "Adjusted.Pvalue",
                                 "name")) %>%
        set_rownames("ensembl_gene_id")
    plot_text$ensembl_gene_id <- NULL

    # When there's time, rework this function internally
    CHBUtils::volcano_density_plot(df,
                                   lfc.cutoff = lfc,
                                   # This isn't documented...
                                   plot_text = plot_text,
                                   pval.cutoff = alpha)

}
