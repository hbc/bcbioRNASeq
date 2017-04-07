#' @rdname qc_plots
#' @description Genes detection saturation plot
#' @return Smooth plot
#' @export
plot_gene_detection_saturation <- function(bcbio, raw_counts) {
    check_bcbio(bcbio)

    plot <- import_summary(bcbio) %>%
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

    show(plot)
}
