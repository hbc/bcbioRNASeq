#' Differential Expression Results Tables
#'
#' @note Log fold change cutoff ("`lfc`") does not apply to statistical
#'   hypothesis testing, only gene filtering in the results tables. See
#'   [DESeq2::results()] for additional information about using `lfcThreshold`
#'   and `altHypothesis` to set an alternative hypothesis based on expected fold
#'   changes.
#'
#' @name resultsTables
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param rowData Join Ensembl gene annotations to the results. Apply gene
#'   identifier to symbol mappings. A previously saved `data.frame` is
#'   recommended. Alternatively if set `NULL`, then gene annotations will not be
#'   added to the results.
#' @param summary Show summary statistics.
#' @param write Write CSV files to disk.
#' @param dir Local directory path where to write the results tables.
#' @param dropboxDir Dropbox directory path where to archive the results tables
#'   for permanent storage (e.g. Stem Cell Commons). When this option is
#'   enabled, unique links per file are generated internally with the rdrop2
#'   package.
#' @param rdsToken RDS file token to use for Dropbox authentication.
#'
#' @return Results `list`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/res_small.rda", package = "bcbioRNASeq"))
#'
#' # DESeqResults ====
#' resTbl <- resultsTables(
#'     object = res_small,
#'     lfc = 0.25,
#'     rowData = rowData(bcb_small),
#'     summary = TRUE,
#'     headerLevel = 2L,
#'     write = FALSE
#' )
#' names(resTbl)
NULL



# Constructors =================================================================
#' Markdown List of Results Files
#'
#' Enables looping of results contrast file links for RMarkdown.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom basejump markdownHeader markdownList
#'
#' @return [writeLines()] output.
#' @noRd
.markdownResultsTables <- function(object, headerLevel = 2L) {
    assert_is_list(object)
    assert_is_subset(
        c("all", "deg", "degLFCDown", "degLFCUp"),
        names(object)
    )
    assertIsImplicitInteger(headerLevel)

    # Prioritze `dropboxFiles` over `localFiles` for path return
    assert_are_intersecting_sets(
        c("dropboxFiles", "localFiles"),
        names(object)
    )
    if ("dropboxFiles" %in% names(object)) {
        paths <- vapply(
            X = object[["dropboxFiles"]],
            FUN = function(x) {
                x[["url"]]
            },
            FUN.VALUE = "character"
        )
        basenames <- basename(paths) %>%
            gsub("\\?.*$", "", .)
    } else if ("localFiles" %in% names(object)) {
        paths <- object[["localFiles"]]
        basenames <- basename(paths)
    }
    names(basenames) <- names(paths)

    markdownHeader("Results tables", level = headerLevel, asis = TRUE)
    markdownList(c(
        paste0(
            "[`", basenames[["all"]], "`]",
            "(", paths[["all"]], "): ",
            "All genes, sorted by Ensembl identifier."
        ),
        paste0(
            "[`", basenames[["deg"]], "`]",
            "(", paths[["deg"]], "): ",
            "Genes that pass the alpha (FDR) cutoff."
        ),
        paste0(
            "[`", basenames[["degLFCUp"]], "`]",
            "(", paths[["degLFCUp"]], "): ",
            "Upregulated DEG; positive log2 fold change."
        ),
        paste0(
            "[`", basenames[["degLFCDown"]], "`]",
            "(", paths[["degLFCDown"]], "): ",
            "Downregulated DEG; negative log2 fold change."
        )
    ), asis = TRUE)
}



#' @importFrom basejump camel initializeDirectory markdownHeader markdownList
#'   sanitizeRowData snake
#' @importFrom bcbioBase copyToDropbox
#' @importFrom dplyr left_join
#' @importFrom readr write_csv
#' @importFrom tibble rownames_to_column
.resultsTables.DESeqResults <- function(  # nolint
    object,
    lfc = 0L,
    rowData = NULL,
    summary = TRUE,
    headerLevel = 2L,
    write = FALSE,
    dir = ".",
    dropboxDir = NULL,
    rdsToken = NULL,
    ...
) {
    # Legacy arguments =========================================================
    call <- match.call(expand.dots = TRUE)
    # annotable
    if ("annotable" %in% names(call)) {
        warn("Use `rowData` instead of `annotable`")
        rowData <- call[["annotable"]]
    }

    # Assert checks ============================================================
    # Passthrough: headerLevel, dropboxDir, rdsToken
    validObject(object)
    assert_is_a_number(lfc)
    assert_all_are_non_negative(lfc)
    assert_is_any_of(rowData, c("data.frame", "NULL"))
    assert_is_a_bool(summary)
    assert_is_a_bool(write)
    dir <- initializeDirectory(dir)

    # Extract internal parameters from DESeqResults object =====================
    contrast <- contrastName(object)
    fileStem <- snake(contrast)
    # Alpha level, slotted in `DESeqResults` metadata
    alpha <- metadata(object)[["alpha"]]
    assert_is_a_number(alpha)

    # Prepare the results tables ===============================================
    all <- object %>%
        as.data.frame() %>%
        rownames_to_column("geneID") %>%
        as("tibble") %>%
        camel(strict = FALSE) %>%
        .[order(.[["geneID"]]), , drop = FALSE]

    # Add Ensembl gene annotations (rowData), if desired
    if (has_dims(rowData)) {
        # Drop the nested lists (e.g. entrezID), otherwise the CSVs will fail
        # to save when `write = TRUE`.
        rowData <- sanitizeRowData(rowData)
        all <- left_join(
            x = all,
            y = as.data.frame(rowData),
            by = "geneID"
        )
    }

    # Check for overall gene expression with base mean
    baseMeanGt0 <- all %>%
        .[order(.[["baseMean"]], decreasing = TRUE), , drop = FALSE] %>%
        .[.[["baseMean"]] > 0L, , drop = FALSE]
    baseMeanGt1 <- baseMeanGt0 %>%
        .[.[["baseMean"]] > 1L, , drop = FALSE]

    # All DEG tables are sorted by BH adjusted P value
    deg <- all %>%
        .[!is.na(.[["padj"]]), , drop = FALSE] %>%
        .[.[["padj"]] < alpha, , drop = FALSE] %>%
        .[order(.[["padj"]]), , drop = FALSE]
    degLFC <- deg %>%
        .[.[["log2FoldChange"]] > lfc |
            .[["log2FoldChange"]] < -lfc, , drop = FALSE]
    degLFCUp <- degLFC %>%
        .[.[["log2FoldChange"]] > 0L, , drop = FALSE]
    degLFCDown <- degLFC %>%
        .[.[["log2FoldChange"]] < 0L, , drop = FALSE]

    list <- list(
        "contrast" = contrast,
        # Cutoffs
        "alpha" = alpha,
        "lfc" = lfc,
        # Tibbles
        "all" = all,
        "deg" = deg,
        "degLFC" = degLFC,
        "degLFCUp" = degLFCUp,
        "degLFCDown" = degLFCDown
    )

    if (isTRUE(summary)) {
        markdownHeader("Summary statistics", level = headerLevel, asis = TRUE)
        markdownList(c(
            paste(nrow(all), "genes in counts matrix"),
            paste("Base mean > 0:", nrow(baseMeanGt0), "genes (non-zero)"),
            paste("Base mean > 1:", nrow(baseMeanGt1), "genes"),
            paste("Alpha cutoff:", alpha),
            paste("LFC cutoff:", lfc, "(applied in tables only)"),
            paste("DEG pass alpha:", nrow(deg), "genes"),
            paste("DEG LFC up:", nrow(degLFCUp), "genes"),
            paste("DEG LFC down:", nrow(degLFCDown), "genes")
        ), asis = TRUE)
    }

    if (isTRUE(write)) {
        tibbles <- c("all", "deg", "degLFCUp", "degLFCDown")

        # Local files (required) ===============================================
        localFiles <- file.path(
            dir,
            paste0(fileStem, "_", snake(tibbles), ".csv.gz")
        )
        names(localFiles) <- tibbles

        # Write the results tibbles to local directory
        invisible(lapply(
            X = seq_along(localFiles),
            FUN = function(a) {
                write_csv(
                    x = get(tibbles[[a]]),
                    path = localFiles[[a]]
                )
            }
        ))

        # Check that writes were successful
        assert_all_are_existing_files(localFiles)

        # Update the list with the file paths
        list[["localFiles"]] <- localFiles

        # Copy to Dropbox (optional) ===========================================
        if (is.character(dropboxDir)) {
            dropboxFiles <- copyToDropbox(
                files = localFiles,
                dir = dropboxDir,
                rdsToken = rdsToken
            )
            assert_is_list(dropboxFiles)
            list[["dropboxFiles"]] <- dropboxFiles
        }

        # Output file information in Markdown format
        .markdownResultsTables(list, headerLevel = headerLevel)
    }

    list
}



# Methods ======================================================================
#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature("DESeqResults"),
    .resultsTables.DESeqResults
)
