#' @rdname qc_plots
#' @description Plot count density
#'
#' @param run \code{bcbio-nextgen} run object
#' @param normalized_counts Normalized counts matrix
#'
#' @return Density plot
#' @export
plot_count_density <- function(
    run,
    normalized_counts) {
    check_run(run)

    name <- deparse(substitute(normalized_counts))
    melted <- melt_log10(run, normalized_counts)

    plot <- ggplot(
        melted,
        aes_(x = ~counts,
             group = ~samplename)
    ) +
        ggtitle(paste("count density:", name)) +
        geom_density() +
        labs(x = expression(log[10]~counts~per~gene),
             y = "density")

    return(plot)
}
