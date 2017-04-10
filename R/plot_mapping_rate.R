#' @rdname qc_plots
#' @description Mapping rate plot
#'
#' @param run \code{bcbio-nextgen} run object
#'
#' @return Bar plot
#' @export
plot_mapping_rate <- function(run) {
    check_run(run)

    plot <- import_summary(run) %>%
        ggplot(
            aes_(x = ~description,
                 y = ~mapped_reads / total_reads * 100,
                 fill = ~qc_color)
        ) +
        ggtitle("mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = "orange",
                   size = 2,
                   yintercept = 70) +
        geom_hline(color = "green",
                   size = 2,
                   yintercept = 90) +
        labs(x = "sample",
             y = "mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()

    return(plot)
}
