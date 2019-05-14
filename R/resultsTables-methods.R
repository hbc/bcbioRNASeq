# Recommend deprecation of this function in favor of using DESeqAnalysis package
# instead in a future update.



#' Differential expression results tables
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
#' @param counts `DESeqDataSet`. Corresponding object containing the counts that
#'   were used to generate the `DESeqResults` object.
#' @param summary `boolean`. Show summary statistics.
#' @param write `boolean`. Write CSV files to disk.
#' @param dropboxDir `string` or `NULL`. Dropbox directory path where to archive
#'   the results tables for permanent storage (e.g. Stem Cell Commons). When
#'   this option is enabled, unique links per file are generated internally with
#'   the rdrop2 package. Note that local files are written to [base::tempdir()]
#'   and the `dir` argument is ignored, if this is enabled.
#' @param rdsToken `string` or `NULL`. RDS file token to use for Dropbox
#'   authentication. If set `NULL` and `dropboxDir` is defined, then an
#'   interactive prompt will appear requesting authorization.
#' @param ... Additional arguments.
#'
#' @return `list` containing modified `DESeqResults` return, including
#'   additional gene-level metadata and normalized counts.
#'
#' @examples
#' # DESeqResults, DESeqDataSet ====
#' x <- resultsTables(results = res_small, counts = dds_small)
#' names(x)
#' x[["deg"]] %>% str()
NULL



#' @rdname resultsTables
#' @export
setGeneric(
    name = "resultsTables",
    def = function(results, counts, ...) {
        standardGeneric("resultsTables")
    }
)



.degList <- function(
    results,
    alpha,
    lfcThreshold = 0L,
    contrast
) {
    assert_is_subset(c("log2FoldChange", "padj"), colnames(results))
    if (missing(alpha)) {
        alpha <- metadata(results)[["alpha"]]
    }
    assert_is_a_number(alpha)
    assert_is_a_number(lfcThreshold)
    if (missing(contrast)) {
        contrast <- contrastName(results)
    }
    assert_is_a_string(contrast)

    all <- as.data.frame(results)
    # DEG tables are sorted by BH adjusted P value
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

    list(
        deg = deg,
        degLFC = degLFC,
        degLFCUp = degLFCUp,
        degLFCDown = degLFCDown,
        all = all,
        contrast = contrast,
        alpha = alpha,
        lfcThreshold = lfcThreshold
    )
}



#' Markdown list of results files
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
        # nocov start : use local Dropbox token
        paths <- vapply(
            X = object[["dropboxFiles"]],
            FUN = function(x) {
                x[["url"]]
            },
            FUN.VALUE = "character"
        )
        basenames <- basename(paths) %>%
            gsub("\\?.*$", "", .)
        # nocov end
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



`resultsTables.DESeqResults,DESeqDataSet` <-  # nolint
    function(
        results,
        counts,
        alpha,
        lfcThreshold = 0L,
        summary = TRUE,
        write = FALSE,
        headerLevel = 2L,
        dir = ".",
        dropboxDir = NULL,
        rdsToken = NULL
    ) {
        validObject(results)
        validObject(counts)
        assert_are_identical(rownames(results), rownames(counts))
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assert_is_a_bool(summary)
        assert_is_a_bool(write)

        # Write local files to tempdir if Dropbox mode is enabled
        if (is_a_string(dropboxDir)) {
            dir <- tempdir()  # nocov
        } else {
            dir <- initDir(dir)
        }

        # Extract internal parameters from DESeqResults object -----------------
        if (missing(alpha)) {
            alpha <- metadata(results)[["alpha"]]
        }
        assert_is_a_number(alpha)
        contrast <- contrastName(results)
        fileStem <- snake(contrast)

        results <- as.data.frame(results)
        results[["geneID"]] <- rownames(results)

        # Add gene annotations (rowData) ---------------------------------------
        rowData <- rowRanges(counts) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            # Drop any nested complex columns (e.g. entrezID list)
            select_if(is.atomic) %>%
            select(!!!syms(setdiff(colnames(.), colnames(results)))) %>%
            column_to_rownames()
        assert_are_identical(rownames(results), rownames(rowData))
        results <- cbind(results, rowData)

        # Add normalized counts matrix -----------------------------------------
        matrix <- counts(counts, normalized = TRUE)
        # Use the `sampleName` metadata for columns, if defined
        if ("sampleName" %in% colnames(colData(counts))) {
            colnames(matrix) <- colData(counts)[["sampleName"]]
            matrix <- matrix[, sort(colnames(matrix)), drop = FALSE]
        }
        assert_are_identical(rownames(results), rownames(matrix))
        results <- cbind(results, matrix)

        # Check for overall gene expression with base mean
        baseMeanGt0 <- results %>%
            .[.[["baseMean"]] > 0L, , drop = FALSE] %>%
            nrow()
        baseMeanGt1 <- results %>%
            .[.[["baseMean"]] > 1L, , drop = FALSE] %>%
            nrow()

        list <- .degList(
            results = results,
            alpha = alpha,
            lfcThreshold = lfcThreshold,
            contrast = contrast
        )

        # Print a markdown header containing the contrast (useful for looping)
        if (isTRUE(summary) || isTRUE(write)) {
            markdownHeader(contrast, level = headerLevel, asis = TRUE)
        }

        if (isTRUE(summary)) {
            markdownHeader(
                "Summary statistics",
                level = headerLevel + 1L,
                asis = TRUE
            )
            markdownList(c(
                paste(nrow(results), "genes in counts matrix"),
                paste("Base mean > 0:", baseMeanGt0, "genes (non-zero)"),
                paste("Base mean > 1:", baseMeanGt1, "genes"),
                paste("Alpha:", alpha),
                paste("LFC threshold:", lfcThreshold),
                paste("DEG pass alpha:", nrow(list[["deg"]]), "genes"),
                paste("DEG LFC up:", nrow(list[["degLFCUp"]]), "genes"),
                paste("DEG LFC down:", nrow(list[["degLFCDown"]]), "genes")
            ), asis = TRUE)
        }

        if (isTRUE(write)) {
            tables <- list[c("all", "deg", "degLFCUp", "degLFCDown")]

            # Local files (required) -------------------------------------------
            localFiles <- file.path(
                dir,
                paste0(fileStem, "_", snake(names(tables)), ".csv.gz")
            )
            names(localFiles) <- names(tables)

            # Write the results tables to local directory
            if (length(dropboxDir)) {
                # nocov start : use local Dropbox token
                message(paste(
                    "Writing",
                    toString(basename(localFiles)),
                    "to Dropbox",
                    paste0("(", dropboxDir, ")")
                ))
                # nocov end
            } else {
                message(paste(
                    "Writing", toString(basename(localFiles)), "to", dir
                ))
            }

            mapply(
                x = tables,
                path = localFiles,
                FUN = function(x, path) {
                    assert_is_subset("geneID", colnames(x))
                    write_csv(x = x, path = path)
                }
            )

            # Check that writes were successful
            assert_all_are_existing_files(localFiles)

            # Update the list with the file paths
            list[["localFiles"]] <- localFiles

            # Copy to Dropbox (optional) ---------------------------------------
            if (is.character(dropboxDir)) {
                # nocov start
                dropboxFiles <- copyToDropbox(
                    files = localFiles,
                    dir = dropboxDir,
                    rdsToken = rdsToken
                )
                assert_is_list(dropboxFiles)
                list[["dropboxFiles"]] <- dropboxFiles
                # nocov end
            }

            # Output file information in Markdown format
            .markdownResultsTables(list, headerLevel = headerLevel + 1L)
        }

        list
    }



#' @rdname resultsTables
#' @export
setMethod(
    f = "resultsTables",
    signature = signature(
        results = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    definition = `resultsTables.DESeqResults,DESeqDataSet`
)



# Minimal method that is used by other functions inside the package.

`resultsTables.DESeqResults,NULL` <-  # nolint
    function(
        results,
        counts = NULL,
        alpha,
        lfcThreshold = 0L
    ) {
        .degList(
            results = results,
            alpha = alpha,
            lfcThreshold = lfcThreshold
        )
    }



#' @rdname resultsTables
#' @export
setMethod(
    f = "resultsTables",
    signature = signature(
        results = "DESeqResults",
        counts = "missingOrNULL"
    ),
    definition = `resultsTables.DESeqResults,NULL`
)
