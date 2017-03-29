#' Genes detection saturation plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param counts Raw counts matrix
#' @param summary \code{bcbio-rnaseq} summary report
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_gene_detection_saturation(summary, raw_counts)
#' }
plot_gene_detection_saturation <- function(summary, counts) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~mapped_reads / 1e6,
                          y = ~colSums(counts > 0),
                          color = ~qc_color,
                          shape = ~qc_color)) +
        ggplot2::ggtitle("gene detection saturation") +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        ggplot2::labs(x = "mapped reads (million)",
                      y = "gene count",
                      color = "",
                      shape = "")
}
