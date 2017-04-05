#' Differentially expressed gene (DEG) tables
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#'
#' @param bcbio bcbio list object
#' @param res \code{DESeqResults} object
#' @param lfc Log fold change ratio (base 2) cutoff. Does not apply to
#'   statistical hypothesis testing, only gene filtering in the results tables.
#'   See \code{?DESeq2::results} for additional information about using
#'   \code{lfcThreshold} and \code{altHypothesis} to set an alternative
#'   hypothesis based on expected fold changes.
#' @param write_csv Write CSV files
#'
#' @return Results list
#' @export
#'
#' @examples
#' \dontrun{
#' deg_tables(bcbio, res, lfc = 0.25)
#' }
deg_tables <- function(
    bcbio,
    res,
    lfc = 0,
    write_csv = TRUE) {
    check_bcbio_object(bcbio)
    if (class(res)[1] != "DESeqResults") {
        stop("DESeqResults required")
    }
    name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res) %>% gsub(" ", "_", .)
    alpha <- res@metadata$alpha

    # Running biomaRt without setting `values = results@rownames` is faster
    annotations <- ensembl_annotations(
        bcbio,
        attributes = c("external_gene_name",
                       "description",
                       "gene_biotype")
    )

    all <- res %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        set_names_snake %>%
        dplyr::left_join(annotations, by = "ensembl_gene_id") %>%
        dplyr::arrange_(.dots = "ensembl_gene_id") %>%
        set_rownames("ensembl_gene_id")

    # All DEG tables sorted by BH adjusted P value
    deg <- all %>%
        dplyr::filter_(.dots = ~padj < alpha) %>%
        dplyr::arrange_(.dots = "padj") %>%
        set_rownames("ensembl_gene_id")
    deg_lfc <- deg %>%
        dplyr::filter_(.dots = ~log2_fold_change < -lfc |
                           log2_fold_change > lfc) %>%
        set_rownames("ensembl_gene_id")
    deg_lfc_up <- deg_lfc %>%
        dplyr::filter_(.dots = ~log2_fold_change > 0) %>%
        set_rownames("ensembl_gene_id")
    deg_lfc_down <- deg_lfc %>%
        dplyr::filter_(.dots = ~log2_fold_change < 0) %>%
        set_rownames("ensembl_gene_id")

    base_mean_gt0 <- dplyr::filter_(all, .dots = ~base_mean > 0)
    base_mean_gt1 <- dplyr::filter_(all, .dots = ~base_mean > 1)

    writeLines(c(
        paste(name, "differential expression tables"),
        paste(nrow(all), "gene annotations"),
        "cutoffs applied",
        paste("    alpha:", alpha),
        paste("    lfc:  ", lfc, "(tables only)"),
        "base mean",
        paste("    > 0:", nrow(base_mean_gt0), "genes (detected)"),
        paste("    > 1:", nrow(base_mean_gt1), "genes"),
        "pass cutoffs",
        paste("    alpha:     ", nrow(deg), "genes"),
        paste("    + lfc up:  ", nrow(deg_lfc_up), "genes"),
        paste("    + lfc down:", nrow(deg_lfc_down), "genes"),
        ""
    ))

    # Write the CSV files to results/
    deg_dir <- file.path("results", "differential_expression")
    dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(write_csv)) {
        write.csv(all,
                  file = file.path(deg_dir,
                                   paste0(contrast_name, "_all_genes.csv")))
        write.csv(deg,
                  file = file.path(deg_dir,
                                   paste0(contrast_name, "_deg.csv")))
        write.csv(deg_lfc_up,
                  file = file.path(deg_dir,
                                   paste0(contrast_name, "_deg_lfc_up.csv")))
        write.csv(deg_lfc_down,
                  file = file.path(deg_dir,
                                   paste0(contrast_name, "_deg_lfc_down.csv")))
    }

    return(list(
        name = name,
        alpha = alpha,
        lfc = lfc,
        all = all,
        deg = deg,
        deg_lfc = deg_lfc,
        deg_lfc_up = deg_lfc_up,
        deg_lfc_down = deg_lfc_down,
        res = res
    ))
}
