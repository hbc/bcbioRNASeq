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
#' @param dds DESeq2 data set
#' @param annotations Annotations data frame containing gene names
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
#' plot_volcano(dds,
#'              annotations = annotations,
#'              alpha = 0.05,
#'              lfc = 1,
#'              text_alpha = 1e-6,
#'              text_lfc = 1)
#' }
plot_volcano <- function(dds,
                         annotations,
                         alpha = 0.05,
                         lfc = 1,
                         text_alpha = 1e-6,
                         text_lfc = 1) {
    df <- dds %>%
        DESeq2::results(.) %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        # Filter zero counts, select columns
        dplyr::filter_(.dots = ~!is.na(log2FoldChange)) %>%
        dplyr::select_(.dots = c("ensembl_gene_id",
                                 "log2FoldChange",
                                 "padj")) %>%
        # Set up colname consistency with `volcano_density_plot()`
        dplyr::rename_(.dots = c(
            "Adjusted.Pvalue" = "padj",
            "logFC" = "log2FoldChange"
        )) %>%
        set_rownames("ensembl_gene_id") %>%
        dplyr::select_(.dots = c("logFC", "Adjusted.Pvalue"))

    plot_text <- df %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        dplyr::left_join(annotations, by = "ensembl_gene_id") %>%
        dplyr::rename_(.dots = c("name" = "external_gene_name")) %>%
        dplyr::arrange_(.dots = "name") %>%
        dplyr::select_(.dots = c("ensembl_gene_id",
                                 "name",
                                 "logFC",
                                 "Adjusted.Pvalue")) %>%
        dplyr::filter_(.dots = ~Adjusted.Pvalue < 1e-6) %>%
        dplyr::filter_(.dots = ~logFC < -1 | logFC > 1) %>%
        set_rownames("ensembl_gene_id") %>%
        # Must supply in this order
        dplyr::select_(.dots = c("logFC", "Adjusted.Pvalue", "name"))

    # When there's time, rework this function internally
    CHBUtils::volcano_density_plot(df,
                                   lfc.cutoff = lfc,
                                   # This isn't documented...
                                   plot_text = plot_text,
                                   pval.cutoff = alpha)

}
