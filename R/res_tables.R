#' Differential expression results tables
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
#' @param res \code{DESeqResults}
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
#' res_tables(run, res, lfc = 0.25)
#' }
res_tables <- function(
    run,
    res,
    lfc = 0) {
    check_run(run)
    check_res(res)

    name <- deparse(substitute(res))
    contrast <- res_contrast_name(res)
    contrast_name <- contrast %>% gsub(" ", "_", .)
    alpha <- res@metadata$alpha
    ensembl <- run$ensembl

    all <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        set_names_snake %>%
        left_join(ensembl, by = "ensembl_gene_id") %>%
        .[order(.$ensembl_gene_id), ]

    # All DEG tables sorted by BH adjusted P value
    deg <- all %>%
        .[.$padj < alpha, ] %>%
        .[order(.$padj), ]
    glimpse(deg)

    deg_lfc <- deg %>%
        .[which(.$log2_fold_change > lfc | .$log2_fold_change < -lfc), ]
    deg_lfc_up <- deg_lfc %>%
        .[.$log2_fold_change > 0, ]
    deg_lfc_down <- deg_lfc %>%
        .[.$log2_fold_change < 0, ]

    base_mean_gt0 <- all %>%
        .[.$base_mean > 0, ] %>%
        .[order(.$base_mean, decreasing = TRUE), ]
    base_mean_gt1 <- base_mean_gt0 %>%
        .[.$base_mean > 1, ]

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

    # Write the CSV files to `results/`
    de_dir <- file.path("results", "de")
    dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

    write_csv(
        all,
        file.path(de_dir,
                  paste(name,
                        contrast_name, "all_genes.csv.gz",
                        sep = "_")))
    write_csv(
        deg,
        file.path(de_dir,
                  paste(name,
                        contrast_name,
                        "deg.csv.gz",
                        sep = "_")))
    write_csv(
        deg_lfc_up,
        file.path(de_dir,
                  paste(name,
                        contrast_name,
                        "deg_lfc_up.csv.gz",
                        sep = "_")))
    write_csv(
        deg_lfc_down,
        file.path(de_dir,
                  paste(name,
                        contrast_name,
                        "deg_lfc_down.csv.gz",
                        sep = "_")))

    return(list(
        name = name,
        contrast = contrast,
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
