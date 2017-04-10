#' @rdname qc_plots
#' @description Genes detection saturation plot
#'
#' @param run \code{bcbio-nextgen} run object
#' @param raw_counts Raw counts matrix. Can be obtained from \code{DESeqDataSet}
#'   by running \code{counts(normalized = FALSE)}.
#'
#' @return Smooth plot
#' @export
plot_gene_detection_saturation <- function(run, raw_counts) {
    check_run(run)

    plot <- import_summary(run) %>%
        ggplot(
            aes_(x = ~mapped_reads / 1e6,
                 y = ~colSums(raw_counts > 0),
                 color = ~qc_color,
                 shape = ~qc_color)) +
        ggtitle("gene detection saturation") +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(x = "mapped reads (million)",
             y = "gene count",
             color = "",
             shape = "")

    return(plot)
}
