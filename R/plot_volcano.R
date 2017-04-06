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
#' @param res \code{DESeqResults} object
#' @param lfc Log fold change ratio (base 2) cutoff for coloring
#' @param text_labels Number of text labels to plot
#'
#' @return Volcano plot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_volcano(bcbio,
#'              res,
#'              alpha = 0.05,
#'              lfc = 1)
#' }
plot_volcano <- function(bcbio,
                         res,
                         lfc = 1,
                         text_labels = 20) {
    check_bcbio(bcbio)
    check_res(res)

    name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)

    alpha <- res@metadata$alpha

    # Prepare data frame for `CHBUtils::volcano_density_plot()`. Note that
    # DESeq2 result table columns must be renamed.
    df <- res %>%
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
        dplyr::arrange_(.dots = "Adjusted.Pvalue") %>%
        set_rownames("ensembl_gene_id") %>%
        dplyr::select_(.dots = c("logFC", "Adjusted.Pvalue"))

    # Automatically label the top genes
    plot_text <- df[1:text_labels, ] %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
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

    # When there's time, rework this function internally in the package
    CHBUtils::volcano_density_plot(
        df,
        lfc.cutoff = lfc,
        # This isn't documented...
        plot_text = plot_text,
        pval.cutoff = alpha,
        title = paste0(name, ": ", contrast_name))
}
