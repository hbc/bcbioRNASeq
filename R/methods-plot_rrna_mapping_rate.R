#' Plot rRNA Mapping Rate
#'
#' @rdname plot_rrna_mapping_rate
#' @author Michael Steinbaugh
#'
#' @export
plot_rrna_mapping_rate <- function(
    bcb,
    warn_limit = 10L,
    interesting_groups = NULL) {
    # FIXME Upgrade to S4 method
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~r_rna_rate * 100L,
             fill = interesting_groups)) +
        ggtitle("rRNA mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "rRNA mapping rate (%)",
             fill = "") +
        coord_flip()
}
