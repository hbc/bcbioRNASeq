#' Differential Expression Results Tables
#'
#' @author Michael Steinbaugh
#'
#' @param bcb [bcbioRNADataSet].
#' @param res [DESeqResults].
#' @param lfc Log fold change ratio (base 2) cutoff. Does not apply to
#'   statistical hypothesis testing, only gene filtering in the results tables.
#'   See [results()] for additional information about using `lfcThreshold` and
#'   `altHypothesis` to set an alternative hypothesis based on expected fold
#'   changes.
#' @param write Write CSV files to disk.
#' @param dir Directory path where to write files.
#'
#' @return Results list.
#' @export
#'
#' @examples
#' \dontrun{
#' resultsTables(bcb, res, lfc = 0.25)
#' }
resultsTables <- function(
    bcb,
    res,
    lfc = 0L,
    write = TRUE,
    dir = file.path("results", "de")) {
    name <- deparse(substitute(res))
    contrast <- .resContrastName(res)
    fileStem <- paste(name,
                       str_replace_all(contrast, " ", "_"),
                       sep = "_")

    # Alpha level, from [DESeqResults]
    alpha <- metadata(res)[["alpha"]]
    annotations <- rowData(bcb) %>%
        as.data.frame %>%
        rename(ensgene = .data[["ensgene"]])

    all <- res %>%
        as.data.frame %>%
        rownames_to_column("ensgene") %>%
        as("tibble") %>%
        camel %>%
        left_join(annotations, by = "ensgene") %>%
        arrange(!!sym("ensgene"))

    # Check for overall gene expression with base mean
    baseMeanGt0 <- all %>%
        arrange(desc(!!sym("baseMean"))) %>%
        .[.[["baseMean"]] > 0L, ]
    baseMeanGt1 <- baseMeanGt0 %>%
        .[.[["baseMean"]] > 1L, ]

    # All DEG tables are sorted by BH adjusted P value
    deg <- all %>%
        .[.[["padj"]] < alpha, ] %>%
        arrange(!!sym("padj"))
    degLFC <- deg %>%
        .[.[["log2FoldChange"]] > lfc |
              .[["log2FoldChange"]] < -lfc, ]
    degLFCUp <- degLFC %>%
        .[.[["log2FoldChange"]] > 0L, ]
    degLFCDown <- degLFC %>%
        .[.[["log2FoldChange"]] < 0L, ]

    resTbl <- list(
        res = res,
        name = name,
        contrast = contrast,
        # Cutoffs
        alpha = alpha,
        lfc = lfc,
        # Tibbles
        all = all,
        deg = deg,
        degLFC = degLFC,
        degLFCUp = degLFCUp,
        degLFCDown = degLFCDown)

    writeLines(c(
        paste(name, "differential expression tables"),
        paste("-", nrow(all), "gene annotations"),
        "- cutoffs applied:",
        paste("    - alpha:", alpha),
        paste("    - lfc:  ", lfc, "(tables only)"),
        "- gene detection:",
        paste("    - base mean > 0:", nrow(baseMeanGt0), "genes"),
        paste("    - base mean > 1:", nrow(baseMeanGt1), "genes"),
        "- pass cutoffs:",
        paste("    - alpha:   ", nrow(deg), "genes"),
        paste("    - lfc up:  ", nrow(degLFCUp), "genes"),
        paste("    - lfc down:", nrow(degLFCDown), "genes")))

    if (isTRUE(write)) {
        # Write the CSV files
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)

        allFile <- paste(fileStem, "all_genes.csv.gz", sep = "_")
        write_csv(all, file.path(dir, allFile))

        degFile <- paste(fileStem, "deg.csv.gz", sep = "_")
        write_csv(deg, file.path(dir, degFile))

        degLFCUpFile <- paste(fileStem, "deg_lfc_up.csv.gz", sep = "_")
        write_csv(degLFCUp, file.path(dir, degLFCUpFile))

        degLFCDownFile <- paste(fileStem, "deg_lfc_down.csv.gz", sep = "_")
        write_csv(degLFCDown, file.path(dir, degLFCDownFile))

        # Output file information in Markdown format
        .mdResultsTables(resTbl, dir)
    }

    resTbl
}
