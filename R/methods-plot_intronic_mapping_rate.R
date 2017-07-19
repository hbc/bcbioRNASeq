#' Plot Intronic Mapping Rate
#'
#' @rdname plot_intronic_mapping_rate
#' @export
plot_intronic_mapping_rate <- function(
    bcb,
    warn_limit = 20L,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~intronic_rate * 100L,
             fill = interesting_groups)) +
        ggtitle("intronic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "intronic mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}
