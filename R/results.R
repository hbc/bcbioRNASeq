#' Print summary statistics of alpha level cutoffs
#'
#' @author Michael Steinbaugh
#'
#' @param dds [DESeqDataSet].
#' @param alpha Numeric vector of desired alpha cutoffs.
#' @param ... Passthrough parameters to [results()].
#'
#' @return Printed summary.
#' @export
alpha_summary <- function(
    dds,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
    ...) {
    name <- deparse(substitute(dds))
    print(name)
    lapply(seq_along(alpha), function(a) {
        writeLines(
            paste(name,
                  paste("alpha", alpha[a], sep = " = "),
                  sep = " : "))
        results(dds, alpha = alpha[a], ...) %>% summary
    }) %>% invisible
}



#' Metadata table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as `intgroup` in the
#' [bcbioRnaDataSet] object.
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioRnaDataSet].
#'
#' @return Data frame containing only the columns of interest.
#' @export
metadata_table <- function(run) {
    run$metadata %>%
        remove_rownames %>%
        kable(caption = "Sample metadata")
}



#' Differential expression results tables
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioRnaDataSet].
#' @param res [DESeqResults].
#' @param lfc Log fold change ratio (base 2) cutoff. Does not apply to
#'   statistical hypothesis testing, only gene filtering in the results tables.
#'   See [results()] for additional information about using `lfcThreshold` and
#'   `altHypothesis` to set an alternative hypothesis based on expected fold
#'   changes.
#' @param write Write CSV files to disk.
#' @param dir Directory path where to write files.
#' @param print Print summary statistics and links to files.
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
    lfc = 0,
    write = TRUE,
    dir = "results/de",
    print = TRUE) {
    check_run(run)
    check_res(res)
    import_tidy_verbs()

    name <- deparse(substitute(res))
    contrast <- res_contrast_name(res)
    file_stem <- paste(name,
                       str_replace_all(contrast, " ", "_"),
                       sep = "_")

    alpha <- res@metadata$alpha
    annotations <- gene_level_annotations(run)

    all <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        snake %>%
        left_join(annotations, by = "ensembl_gene_id") %>%
        arrange(!!sym("ensembl_gene_id"))

    # Check for overall gene expression with base mean
    base_mean_gt0 <- all %>%
        arrange(desc(!!sym("base_mean"))) %>%
        filter(.data$base_mean > 0)
    base_mean_gt1 <- base_mean_gt0 %>%
        filter(.data$base_mean > 1)

    # All DEG tables are sorted by BH adjusted P value
    deg <- all %>%
        filter(.data$padj < alpha) %>%
        arrange(!!sym("padj"))
    deg_lfc <- deg %>%
        filter(.data$log2_fold_change > lfc | .data$log2_fold_change < -lfc)
    deg_lfc_up <- deg_lfc %>%
        filter(.data$log2_fold_change > 0)
    deg_lfc_down <- deg_lfc %>%
        filter(.data$log2_fold_change < 0)

    # Write the CSV files to `results/`
    if (isTRUE(write)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)

        all_file <- paste(file_stem, "all_genes.csv.gz", sep = "_")
        write_csv(all, file.path(dir, all_file))

        deg_file <- paste(file_stem, "deg.csv.gz", sep = "_")
        write_csv(deg, file.path(dir, deg_file))

        deg_lfc_up_file <- paste(file_stem, "deg_lfc_up.csv.gz", sep = "_")
        write_csv(deg_lfc_up, file.path(dir, deg_lfc_up_file))

        deg_lfc_down_file <- paste(file_stem, "deg_lfc_down.csv.gz", sep = "_")
        write_csv(deg_lfc_down, file.path(dir, deg_lfc_down_file))
    }

    res_tbl <- list(
        res = res,
        name = name,
        contrast = contrast,
        # Cutoffs
        alpha = alpha,
        lfc = lfc,
        # Tibbles
        all = all,
        deg = deg,
        deg_lfc = deg_lfc,
        deg_lfc_up = deg_lfc_up,
        deg_lfc_down = deg_lfc_down,
        # File output
        dir = dir,
        file_stem = file_stem,
        all_file = all_file,
        deg_file = deg_file,
        deg_lfc_up_file = deg_lfc_up_file,
        deg_lfc_down_file = deg_lfc_down_file)

    # Print summary statistics and file paths
    if (isTRUE(print)) {
        writeLines(c(
            paste(name, "differential expression tables"),
            "",
            paste("-", nrow(all), "gene annotations"),
            "- cutoffs applied:",
            paste("    - alpha:", alpha),
            paste("    - lfc:  ", lfc, "(tables only)"),
            "- gene detection:",
            paste("    - base mean > 0:", nrow(base_mean_gt0), "genes"),
            paste("    - base mean > 1:", nrow(base_mean_gt1), "genes"),
            "- pass cutoffs:",
            paste("    - alpha:   ", nrow(deg), "genes"),
            paste("    - lfc up:  ", nrow(deg_lfc_up), "genes"),
            paste("    - lfc down:", nrow(deg_lfc_down), "genes"),
            "",
            ""
        ))
        md_res_tables(res_tbl)
    }

    res_tbl
}



#' Top tables of differential expression results
#'
#' @author Michael Steinbaugh
#'
#' @param res_tbl Results tables generated by [res_tables()].
#' @param n Number genes to report.
#' @param coding Whether to only return coding genes.
#'
#' @return Top tables list.
#' @export
top_tables <- function(
    res_tbl,
    n = 50,
    coding = FALSE) {
    import_tidy_verbs()
    subset_top <- function(df) {
        # Filter for coding genes only, if desired
        if (isTRUE(coding)) {
            df <- filter(df, .data$broad_class == "coding")
        }

        df <- df %>%
            # [top_n()] defaults to last column. `wt` overrides. Here we want
            # to rank by BH adjusted P value (`padj`).
            top_n(n = n, wt = !!sym("padj")) %>%
            mutate(
                base_mean = round(.data$base_mean),
                log2_fold_change = format(.data$log2_fold_change, digits = 3),
                padj = format(.data$padj, digits = 3, scientific = TRUE)) %>%
            select(!!!syms(c("ensembl_gene_id",
                             "base_mean",
                             "log2_fold_change",
                             "padj",
                             "external_gene_name",
                             "broad_class"))) %>%
            remove_rownames

        df
    }

    up <- subset_top(res_tbl$deg_lfc_up)
    down <- subset_top(res_tbl$deg_lfc_down)

    # Captions
    name <- res_tbl$name
    contrast <- res_tbl$contrast
    sep <- " : "
    name_prefix <- paste(name, contrast, sep = sep)

    show(kable(up, caption = paste(name_prefix, "upregulated", sep = sep)))
    show(kable(down, caption = paste(name_prefix, "downregulated", sep = sep)))
}
