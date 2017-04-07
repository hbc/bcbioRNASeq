#' @rdname qc_plots
#' @description Exonic mapping rate plot
#' @return Bar plot
#' @export
plot_exonic_mapping_rate <- function(bcbio) {
    check_bcbio(bcbio)

    plot <- import_summary(bcbio) %>%
        ggplot(
            aes_(x = ~description,
                 # Multiple by 100 here for percentage
                 y = ~exonic_rate * 100,
                 fill = ~qc_color)
        ) +
        ggtitle("exonic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = "green",
                   size = 2,
                   yintercept = 60) +
        labs(x = "sample",
             y = "exonic mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()

    show(plot)
}
