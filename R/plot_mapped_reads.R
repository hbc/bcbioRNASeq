#' @rdname qc_plots
#' @description Mapped reads plot
#' @return Bar plot
#' @export
plot_mapped_reads <- function(bcbio) {
    check_bcbio(bcbio)

    plot <- import_summary(bcbio) %>%
        ggplot(
            aes_(x = ~description,
                 y = ~mapped_reads / 1e6,
                 fill = ~qc_color)
        ) +
        ggtitle("mapped reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = "orange",
                   size = 2,
                   yintercept = 10) +
        geom_hline(color = "green",
                   size = 2,
                   yintercept = 20) +
        labs(x = "sample",
             y = "mapped reads (million)",
             fill = "") +
        coord_flip()

    show(plot)
}
