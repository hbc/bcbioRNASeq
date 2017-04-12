#' Differential expression plots
#'
#' @rdname de_plots
#' @author Michael Steinbaugh



#' @rdname de_plots
#' @description Wrapper for \code{DESeq2::plotMA} that generates a title
#'   automatically
#'
#' @import DESeq2
#'
#' @param res \code{DESeqResults}
#' @param ylim Y-axis maximum (single integer)
#'
#' @return MA plot
#' @export
plot_ma <- function(res, ylim = 2) {
    check_res(res)

    name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)

    plotMA(
        res,
        main = paste(name, contrast_name, sep = " : "),
        ylim = c(-ylim, ylim))
}



#' @rdname de_plots
#' @description Volcano plot
#'
#' @import DESeq2
#' @import dplyr
#' @importFrom CHBUtils volcano_density_plot
#' @importFrom tibble rownames_to_column
#'
#' @param run \code{bcbio-nextgen} run
#' @param lfc Log fold change ratio (base 2) cutoff for coloring
#' @param text_labels Number of text labels to plot
#'
#' @return Volcano plot
#' @export
plot_volcano <- function(
    run,
    res,
    lfc = 1,
    text_labels = 30) {
    check_run(run)
    check_res(res)

    name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)

    alpha <- res@metadata$alpha

    # Prepare data frame for `CHBUtils::volcano_density_plot`. Note that
    # DESeq2 result table columns must be renamed.
    df <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        # Filter zero counts, select columns
        filter_(.dots = ~!is.na(log2FoldChange)) %>%
        select_(.dots = c("ensembl_gene_id",
                          "log2FoldChange",
                          "padj")) %>%
        rename_(.dots = c(
            "Adjusted.Pvalue" = "padj",
            "logFC" = "log2FoldChange"
        )) %>%
        arrange_(.dots = "Adjusted.Pvalue") %>%
        set_rownames(.$ensembl_gene_id) %>%
        select_(.dots = c("logFC", "Adjusted.Pvalue"))

    # Automatically label the top genes
    plot_text <- df[1:text_labels, ] %>%
        rownames_to_column("ensembl_gene_id") %>%
        left_join(
            ensembl_annotations(run, values = .$ensembl_gene_id),
            by = "ensembl_gene_id"
        ) %>%
        rename_(.dots = c("name" = "external_gene_name")) %>%
        select_(.dots = c("ensembl_gene_id",
                          "logFC",
                          "Adjusted.Pvalue",
                          "name")) %>%
        set_rownames(.$ensembl_gene_id)
    plot_text$ensembl_gene_id <- NULL

    # When there's time, rework this function internally in the package
    # https://github.com/hbc/CHBUtils/blob/master/R/volcanoPlot.R
    # `volcano_density_plot()` requires a `data.frame` with two columns:
    # `logFC` and `Adjusted.Pvalue`
    volcano_density_plot(
        df,
        lfc.cutoff = lfc,
        # This isn't documented...
        plot_text = plot_text,
        pval.cutoff = alpha,
        title = paste(name, contrast_name, sep = " : "))
}
