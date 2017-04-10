#' @rdname qc_plots
#' @description rRNA contamination mapping rate
#'
#' @param run \code{bcbio-nextgen} run object
#'
#' @return Bar plot
#' @export
plot_rrna_mapping_rate <- function(run) {
    check_run(run)

    plot <- import_summary(run) %>%
        ggplot(
            aes_(x = ~description,
                 y = ~rrna_rate * 100,
                 fill = ~qc_color)) +
        ggtitle("rRNA mapping rate") +
        geom_bar(stat = "identity") +
        #` geom_hline(linetype = 2, yintercept = 5) +
        geom_hline(color = "orange",
                   size = 2,
                   yintercept = 10) +
        labs(x = "sample",
             y = "rRNA mapping rate (%)",
             fill = "") +
        coord_flip()

    return(plot)
}
