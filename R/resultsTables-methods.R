# FIXME Define and recommend DESeqAnalysis method
# FIXME Add as many assert checks as possible here to look for mismatch.
# FIXME Require valid names for the rownames.
# FIXME Improve the documentation about what types of counts to use here.
# FIXME Recommend using normalized counts or DESeqTransform.
# FIXME Define resultsTables S4 class
# FIXME Simplify how we're handling rowData?



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
#'
#' @return `DESeqResultsTables`.
#'
#' @examples
#' # DESeqResults ====
#' x <- resultsTables(res_small)
#' class(x)
#' slotNames(x)
NULL



# FIXME We need to tweak this function.

#' Markdown List of Results Files
#'
#' Enables looping of results contrast file links for RMarkdown.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @return [writeLines()] output.
#' @noRd
.markdownTables <- function(object, headerLevel = 2L) {
    assert_is_list(object)
    assert_is_subset(
        c("all", "deg", "degDown", "degUp"),
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
            "Genes that pass the alpha (FDR) and",
            "log2 fold change (LFC) cutoffs."
        ),
        paste0(
            "[`", basenames[["degUp"]], "`]",
            "(", paths[["degLFCUp"]], "): ",
            "Upregulated DEG; positive fold change."
        ),
        paste0(
            "[`", basenames[["degDown"]], "`]",
            "(", paths[["degDown"]], "): ",
            "Downregulated DEG; negative fold change."
        )
    ), asis = TRUE)
}



#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature("DESeqResults"),
    function(object) {
        validObject(object)
        assert_is_all_of(object, "DESeqResults")
        assert_is_subset(c("log2FoldChange", "padj"), colnames(object))
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)

        # Set LFC and test (P value) columns.
        lfcCol <- "log2FoldChange"
        testCol <- "padj"
        lfc <- sym(lfcCol)
        test <- sym(testCol)
        assert_is_subset(
            x = c(lfcCol, testCol),
            y = colnames(object)
        )

        # DEG tables are sorted by adjusted P value.
        deg <- object %>%
            as("tbl_df") %>%
            # Remove genes without an adjusted P value.
            filter(!is.na(!!test)) %>%
            # Remove genes that don't pass alpha cutoff.
            filter(!!test < !!alpha) %>%
            # Arrange by adjusted P value.
            arrange(!!test) %>%
            # Remove genes that don't pass LFC threshold.
            filter(!!lfc > !!lfcThreshold | !!lfc < -UQ(lfcThreshold))
        # Get directional subsets.
        degUp <- filter(deg, !!lfc > 0L)
        degDown <- filter(deg, !!lfc < 0L)

        new(
            Class = "DESeqResultsTables",
            all = object,
            deg = as(deg, "DataFrame"),
            degUp = as(degUp, "DataFrame"),
            degDown = as(degDown, "DataFrame")
        )
    }
)



# TODO Using the normalized counts, not the DESeqTransform here.

#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature("DESeqAnalysis"),
    function(
        object,
        summary = TRUE,
        write = FALSE,
        headerLevel = 2L,
        dir = ".",
        dropboxDir = NULL,
        rdsToken = NULL
    ) {
        validObject(object)
        assert_are_identical(rownames(object), rownames(counts))
        assert_is_a_bool(summary)
        assert_is_a_bool(write)

        # Write local files to tempdir if Dropbox mode is enabled
        if (is_a_string(dropboxDir)) {
            dir <- tempdir()  # nocov
        } else {
            dir <- initializeDirectory(dir)
        }

        # Extract internal parameters from DESeqResults object -----------------
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        contrast <- contrastName(object)
        fileStem <- snake(contrast)

        # FIXME Rework this step...
        results <- as.data.frame(object)
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

        # FIXME
        list <- .resultsTables(
            results = results,
            alpha = alpha,
            lfcThreshold = lfcThreshold
        )

        # Print a markdown header containing the contrast (useful for looping).
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

            # Check that writes were successful.
            assert_all_are_existing_files(localFiles)

            # Update the list with the file paths.
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

            # Output file information in Markdown format.
            .markdownTables(list, headerLevel = headerLevel + 1L)
        }

        list
    }
)
