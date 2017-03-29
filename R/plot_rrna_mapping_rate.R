#' rRNA contamination mapping rate
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
#' plot_rrna_mapping_rate(bcbio)
#' }
plot_rrna_mapping_rate <- function(bcbio) {
    import_summary(bcbio) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = ~rrna_rate * 100,
                          fill = ~qc_color)) +
        ggplot2::ggtitle("rRNA mapping rate") +
        ggplot2::geom_bar(stat = "identity") +
        #` ggplot2::geom_hline(linetype = 2, yintercept = 5) +
        ggplot2::geom_hline(color = "red",
                            size = 2,
                            yintercept = 10) +
        ggplot2::labs(x = "sample",
                      y = "rRNA mapping rate (%)",
                      fill = "") +
        ggplot2::coord_flip()
}
