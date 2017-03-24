#' Genes detection saturation plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param raw_counts Raw counts matrix
#' @param summary \code{bcbio-rnaseq} summary report
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_gene_detection_saturation(raw_counts, summary)
#' }
plot_gene_detection_saturation <- function(raw_counts, summary) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~mapped_reads / 1e6,
                          y = ~colSums(raw_counts > 0),
                          color = ~group,
                          shape = ~group)) +
        ggplot2::ggtitle("gene detection saturation") +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        #` expand_limits(x = 0, y = 0) +
        ggplot2::xlab("mapped reads (million)") +
        ggplot2::ylab("gene count")
}
