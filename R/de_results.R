#' Differential gene expression results
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#' @importFrom basejump setNamesSnake
#'
#' @param bcbio bcbio list object
#' @param dds \code{DESeqDataSet} object
#' @param alpha Alpha level cutoff. Used for Wald test BH-adjusted P values.
#' @param lfc Log fold change ratio (base 2) cutoff. Does not apply to
#'   statistical hypothesis testing, only gene filtering in the results tables.
#'   See \code{?DESeq2::results} for additional information about using
#'   \code{lfcThreshold} and \code{altHypothesis} to set an alternative
#'   hypothesis based on expected fold changes.
#'
#' @return Results list
#' @export
#'
#' @examples
#' \dontrun{
#' de_results(bcbio, dds, alpha = 0.1, lfc = 0.25)
#' }
de_results <- function(
    bcbio,
    dds,
    alpha = 0.1,
    lfc = 0) {
    check_bcbio_object(bcbio)
    if (class(dds)[1] != "DESeqDataSet") {
        stop("DESeqDataSet required")
    }
    name <- deparse(substitute(dds))

    # Running biomaRt without setting `values = results@rownames` is faster
    annotations <- ensembl_annotations(
        bcbio,
        attributes = c("external_gene_name",
                       "description",
                       "gene_biotype")
    )

    # Nonzero genes
    counts <- DESeq2::counts(dds, normalized = TRUE)
    nonzero_counts <- counts[rowSums(counts) > 0, ]
    ten_counts <- counts[rowSums(counts) > 10, ]

    # DESeqResults object
    res <- DESeq2::results(dds, alpha = alpha)

    message(name)
    message(paste0("alpha = ", res@metadata$alpha))
    message(paste0("lfc = ", lfc, " (only applied to tables)"))

    all <- res %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        tibble::as_tibble(.) %>%
        basejump::setNamesSnake(.) %>%
        dplyr::left_join(annotations, by = "ensembl_gene_id") %>%
        dplyr::arrange_(.dots = "padj")

    deg <- all %>%
        dplyr::filter_(.dots = ~padj < alpha)

    deg_lfc <- deg %>%
        dplyr::filter_(.dots = ~log2_fold_change < -lfc |
                           log2_fold_change > lfc)

    deg_lfc_up <- deg_lfc %>%
        dplyr::filter_(.dots = ~log2_fold_change > 0) %>%
        # ~-log_fold_change also works.
        # Adding `quote()` here is more readable.
        dplyr::arrange_(.dots = quote(-log2_fold_change))

    deg_lfc_down <- deg_lfc %>%
        dplyr::filter_(.dots = ~log2_fold_change < 0) %>%
        dplyr::arrange_(.dots = "log2_fold_change")

    writeLines(c(
        paste(name, "differential expression"),
        paste(nrow(all), "gene annotations"),
        "cutoffs applied",
        paste("  alpha:", alpha),
        paste("  lfc:  ", lfc, "(tables only)"),
        "counts",
        paste("  >= 1: ", nrow(nonzero_counts), "genes"),
        paste("  >= 10:", nrow(ten_counts), "genes"),
        "significance",
        paste("  alpha:", nrow(deg), "genes"),
        paste("  up:   ", nrow(deg_lfc_up), "genes (lfc applied)"),
        paste("  down: ", nrow(deg_lfc_down), "genes (lfc applied)")
    ))

    return(list(
        name = name,
        alpha = alpha,
        lfc = lfc,
        dds = dds,
        res = res,
        all = all,
        deg = deg,
        deg_lfc = deg_lfc,
        deg_lfc_up = deg_up,
        deg_lfc_down = deg_down
    ))
}
