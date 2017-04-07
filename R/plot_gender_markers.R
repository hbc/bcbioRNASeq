#' @rdname qc_plots
#' @description Plot sexually dimorphic gender marker genes
#'
#' @import dplyr
#' @import ggplot2
#' @import readr
#' @import reshape2
#' @import tibble
#'
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   \code{DESeqDataSet} by running \code{counts(normalized = TRUE)}.
#'   Transcripts per million (TPM) are also acceptable.
#'
#' @return Scatter plot
#' @export
plot_gender_markers <- function(bcbio, normalized_counts) {
    check_bcbio(bcbio)

    name <- deparse(substitute(normalized_counts))
    organism <- bcbio$organism

    # Download the CSV file from seqcloud repo
    csv <- file.path("https://raw.githubusercontent.com",
                     "steinbaugh",
                     "seqcloud",
                     "master",
                     "workflows",
                     "illumina_rnaseq",
                     organism,
                     "gender_markers.csv") %>%
        read_csv(col_types = cols(),
                 na = c("", "#N/A")) %>%
        .[.$include == TRUE, ] %>%
        arrange_(.dots = c("chromosome", "gene_symbol"))

    # Ensembl identifiers
    identifier <- csv[, "ensembl_gene"][[1]] %>% sort %>% unique

    # Scatterplot
    plot <- normalized_counts[identifier, ] %>%
        as.data.frame %>%
        rownames_to_column %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1) %>%
        set_names(c("ensembl_gene",
                    "description",
                    "counts")) %>%
        left_join(csv, by = "ensembl_gene") %>%
        ggplot(
            aes_(~gene_symbol,
                 ~counts,
                 color = ~description,
                 shape = ~chromosome)
        ) +
        ggtitle("gender markers") +
        geom_jitter(size = 4) +
        expand_limits(y = 0) +
        labs(x = "gene",
             y = name) +
        theme(legend.position = "none")

    show(plot)
}
