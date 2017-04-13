#' Differential expression plots
#'
#' @rdname de_plots
#' @author Michael Steinbaugh



#' @rdname de_plots
#' @description Wrapper for \code{DESeq2::plotMA} that generates a title
#'   automatically
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

    # Prepare data frame for `volcano_density_plot`
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

    volcano_density_plot(
        df,
        lfc_cutoff = lfc,
        plot_text = plot_text,
        pval_cutoff = alpha,
        title = paste(name, contrast_name, sep = " : "))
}
