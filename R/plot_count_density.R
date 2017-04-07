#' @rdname qc_plots
#' @description Plot count density
#' @return Density plot
#' @export
plot_count_density <- function(
    bcbio,
    counts) {
    check_bcbio(bcbio)
    name <- deparse(substitute(counts))

    melted <- melt_log10(bcbio, counts)

    plot <- ggplot(
        melted,
        aes_(x = ~counts,
             group = ~samplename)
    ) +
        ggtitle(paste("count density:", name)) +
        geom_density() +
        labs(x = expression(log[10]~counts~per~gene),
             y = "density")

    show(plot)
}
