#' Genes detected plot
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
#' plot_genes_detected(raw_counts, summary)
#' }
plot_genes_detected <- function(raw_counts, summary) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = colSums(raw_counts > 0),
                          fill = ~intgroup)
        ) +
        ggplot2::ggtitle("genes detected") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 20000) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("gene count") +
        ggplot2::coord_flip()
}
