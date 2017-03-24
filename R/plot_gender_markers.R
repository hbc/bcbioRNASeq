#' Graph sexually dimorphic marker genes
#'
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import ggplot2
#' @import readr
#' @import reshape2
#' @import tibble
#'
#' @param counts Counts matrix. Transcripts per million (TPM) or normalized
#'   counts are preferred.
#' @param organism Organism. \code{hsapiens} and \code{mmusculus} are supported.
#'
#' @return Scatterplot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_gender_markers(tpm, "mmusculus")
#' }
plot_gender_markers <- function(counts, organism) {
    # Download the CSV file from seqcloud repo
    csv <- file.path("https://raw.githubusercontent.com",
                     "steinbaugh",
                     "seqcloud",
                     "master",
                     "workflows",
                     "illumina_rnaseq",
                     organism,
                     "gender_markers.csv") %>%
        readr::read_csv(.,
                        col_types = readr::cols(),
                        na = c("", "#N/A")) %>%
        .[.$include == TRUE, ] %>%
        dplyr::arrange_(.dots = c("chromosome",
                                  "gene_symbol"))
    # Ensembl identifiers
    identifier <- csv[, "ensembl_gene"][[1]] %>% sort %>% unique
    # Scatterplot
    counts[identifier, ] %>%
        as.data.frame %>%
        tibble::rownames_to_column(.) %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        reshape2::melt(., id = 1) %>%
        set_names(c("ensembl_gene",
                   "description",
                   "counts")) %>%
        dplyr::left_join(csv, by = "ensembl_gene") %>%
        ggplot2::ggplot(
            ggplot2::aes_(~gene_symbol,
                          ~counts,
                          color = ~description,
                          shape = ~chromosome)
        ) +
        ggplot2::ggtitle("gender markers") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_jitter(size = 4) +
        ggplot2::expand_limits(y = 0) +
        ggplot2::xlab("gene") +
        ggplot2::ylab("tpm")
}
