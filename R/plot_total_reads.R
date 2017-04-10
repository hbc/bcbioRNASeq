#' Quality control plots
#'
#' @author Michael Steinbaugh
#'
#' @rdname qc_plots
#' @description Plot total reads
#'
#' @import ggplot2
#'
#' @param run \code{bcbio-nextgen} run object
#'
#' @return Bar plot
#' @export
plot_total_reads <- function(run) {
    check_run(run)

    plot <- import_summary(run) %>%
        ggplot(aes_(x = ~description,
                    y = ~total_reads / 1e6,
                    fill = ~qc_color)
        ) +
        ggtitle("total reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = "orange",
                   size = 2,
                   yintercept = 10) +
        geom_hline(color = "green",
                   size = 2,
                   yintercept = 20) +
        labs(x = "sample",
             y = "total reads (million)",
             fill = "") +
        coord_flip()

    return(plot)
}
