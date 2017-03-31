#' Mapping rate plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param bcbio bcbio list object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_mapping_rate(bcbio)
#' }
plot_mapping_rate <- function(bcbio) {
    check_bcbio_object(bcbio)
    import_summary(bcbio) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = ~mapped_reads / total_reads * 100,
                          fill = ~qc_color)
        ) +
        ggplot2::ggtitle("mapping rate") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 70) +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 90) +
        ggplot2::labs(x = "sample",
                      y = "mapping rate (%)",
                      fill = "") +
        ggplot2::ylim(0, 100) +
        ggplot2::coord_flip()
}
