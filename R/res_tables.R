#' Differential expression results tables
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param res \linkS4class{DESeqResults}.
#' @param lfc Log fold change ratio (base 2) cutoff. Does not apply to
#'   statistical hypothesis testing, only gene filtering in the results tables.
#'   See [DESeq2::results()] for additional information about using
#'   `lfcThreshold` and `altHypothesis` to set an alternative hypothesis based
#'   on expected fold changes.
#'
#' @return Results list.
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
    import_tidy_verbs()

    name <- deparse(substitute(res))
    contrast <- res_contrast_name(res)
    contrast_name <- contrast %>% gsub(" ", "_", .)
    alpha <- res@metadata$alpha
    annotations <- gene_level_annotations(run)

    all <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        snake %>%
        left_join(annotations, by = "ensembl_gene_id") %>%
        arrange(!!sym("ensembl_gene_id"))

    # All DEG tables sorted by BH adjusted P value
    deg <- all %>%
        filter(.data$padj < alpha) %>%
        arrange(!!sym("padj"))
    glimpse(deg)

    deg_lfc <- deg %>%
        .[which(.$log2_fold_change > lfc | .$log2_fold_change < -lfc), ]
    deg_lfc_up <- deg_lfc %>%
        filter(.data$log2_fold_change > 0)
    deg_lfc_down <- deg_lfc %>%
        filter(.data$log2_fold_change < 0)

    base_mean_gt0 <- all %>%
        arrange(desc(!!sym("base_mean")))
        filter(.data$base_mean > 0)
    base_mean_gt1 <- base_mean_gt0 %>%
        filter(.data$base_mean > 1)

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

    list(name = name,
         contrast = contrast,
         alpha = alpha,
         lfc = lfc,
         all = all,
         deg = deg,
         deg_lfc = deg_lfc,
         deg_lfc_up = deg_lfc_up,
         deg_lfc_down = deg_lfc_down,
         res = res)
}
