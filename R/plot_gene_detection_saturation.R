#' Genes detection saturation plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param bcbio bcbio list object
#' @param counts Raw counts matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_gene_detection_saturation(bcbio, raw_counts)
#' }
plot_gene_detection_saturation <- function(bcbio, counts) {
    check_bcbio(bcbio)
    plot <- import_summary(bcbio) %>%
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
    print(plot)
}
