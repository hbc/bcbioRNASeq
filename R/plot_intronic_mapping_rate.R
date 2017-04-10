#' @rdname qc_plots
#' @description Intronic mapping rate plot
#'
#' @param run \code{bcbio-nextgen} run object
#'
#' @return Bar plot
#' @export
plot_intronic_mapping_rate <- function(run) {
    check_run(run)

    plot <- import_summary(run) %>%
        ggplot(
            aes_(x = ~description,
                 # Multiple by 100 here for percentage
                 y = ~intronic_rate * 100,
                 fill = ~qc_color)
        ) +
        ggtitle("intronic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = "orange",
                   size = 2,
                   yintercept = 20) +
        labs(x = "sample",
             y = "intronic mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()

    return(plot)
}
