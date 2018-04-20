#' Differential Expression Results Tables
#'
#' @note Log fold change cutoff threshold ("`lfcThreshold`") does not apply to
#'   statistical hypothesis testing, only gene filtering in the results tables.
#'   See [DESeq2::results()] for additional information about using
#'   `lfcThreshold` and `altHypothesis` to set an alternative hypothesis based
#'   on expected fold changes.
#'
#' @name resultsTables
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param counts `DESeqDataSet` containing the counts that were used to generate
#'   the `DESeqResults`.
#' @param summary Show summary statistics.
#' @param write Write CSV files to disk.
#' @param dropboxDir Dropbox directory path where to archive the results tables
#'   for permanent storage (e.g. Stem Cell Commons). When this option is
#'   enabled, unique links per file are generated internally with the rdrop2
#'   package.
#' @param rdsToken RDS file token to use for Dropbox authentication.
#'
#' @return `list` containing modified `DESeqResults` return, including
#'   additional gene-level metadata and normalized counts.
#'
#' @examples
#' # DESeqResults, DESeqDataSet ====
#' x <- resultsTables(
#'     results = res_small,
#'     counts = dds_small,
#'     lfcThreshold = 0.25,
#'     summary = TRUE,
#'     write = FALSE
#' )
#' names(x)
NULL



# Constructors =================================================================
#' Markdown List of Results Files
#'
#' Enables looping of results contrast file links for RMarkdown.
#'
#' @author Michael Steinbaugh
#' @keywords internal
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



# Methods ======================================================================
#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature(
        results = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    function(
        results,
        counts,
        lfcThreshold = 0L,
        summary = TRUE,
        write = FALSE,
        headerLevel = 2L,
        dir = ".",
        dropboxDir = NULL,
        rdsToken = NULL
    ) {
        # Passthrough: headerLevel, dropboxDir, rdsToken
        validObject(results)
        validObject(counts)
        assert_are_identical(rownames(results), rownames(counts))
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assert_is_a_bool(summary)
        assert_is_a_bool(write)
        dir <- initializeDirectory(dir)

        # Extract internal parameters from DESeqResults object =================
        contrast <- contrastName(results)
        fileStem <- snake(contrast)
        # Alpha level, slotted in `DESeqResults` metadata
        alpha <- metadata(results)[["alpha"]]
        assert_is_a_number(alpha)

        # Prepare the results tables ===========================================
        all <- results %>%
            as.data.frame() %>%
            camel()

        # Add gene annotations (rowData)
        rowRanges <- rowRanges(counts)
        mcols <- mcols(rowRanges)
        mcols <- mcols[, vapply(
            X = mcols,
            FUN = function(x) {
                is.character(x) || is.factor(x)
            },
            FUN.VALUE = logical(1L)
        )]
        mcols(rowRanges) <- mcols
        rowData <- as.data.frame(rowRanges)
        # Drop columns already present in results
        rowData <- rowData[, setdiff(colnames(rowData), colnames(results))]

        # Now safe to bind the rowData annotations
        all <- cbind(all, rowData)

        # Add the DESeq2 normalized counts
        all <- cbind(all, counts(counts, normalized = TRUE))

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
            .[.[["log2FoldChange"]] > lfcThreshold |
                  .[["log2FoldChange"]] < -lfcThreshold, , drop = FALSE]
        degLFCUp <- degLFC %>%
            .[.[["log2FoldChange"]] > 0L, , drop = FALSE]
        degLFCDown <- degLFC %>%
            .[.[["log2FoldChange"]] < 0L, , drop = FALSE]

        list <- list(
            "contrast" = contrast,
            # Cutoffs
            "alpha" = alpha,
            "lfcThreshold" = lfcThreshold,
            # Data frames
            "all" = all,
            "deg" = deg,
            "degLFC" = degLFC,
            "degLFCUp" = degLFCUp,
            "degLFCDown" = degLFCDown
        )

        if (isTRUE(summary)) {
            markdownHeader(
                "Summary statistics",
                level = headerLevel,
                asis = TRUE
            )
            markdownList(c(
                paste(nrow(all), "genes in counts matrix"),
                paste("Base mean > 0:", nrow(baseMeanGt0), "genes (non-zero)"),
                paste("Base mean > 1:", nrow(baseMeanGt1), "genes"),
                paste("Alpha cutoff:", alpha),
                paste("LFC cutoff:", lfcThreshold, "(applied in tables only)"),
                paste("DEG pass alpha:", nrow(deg), "genes"),
                paste("DEG LFC up:", nrow(degLFCUp), "genes"),
                paste("DEG LFC down:", nrow(degLFCDown), "genes")
            ), asis = TRUE)
        }

        if (isTRUE(write)) {
            tables <- c("all", "deg", "degLFCUp", "degLFCDown")

            # Local files (required) ===========================================
            localFiles <- file.path(
                dir,
                paste0(fileStem, "_", snake(tables), ".csv.gz")
            )
            names(localFiles) <- tables

            # Write the results tables to local directory
            invisible(lapply(
                X = seq_along(localFiles),
                FUN = function(a) {
                    write_csv(
                        x = get(tables[[a]]),
                        path = localFiles[[a]]
                    )
                }
            ))

            # Check that writes were successful
            assert_all_are_existing_files(localFiles)

            # Update the list with the file paths
            list[["localFiles"]] <- localFiles

            # Copy to Dropbox (optional) =======================================
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
)
