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
#' @export
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
#' res <- deseq_small@lfcShrink[[1L]]
#' x <- resultsTables(res)
#' class(x)
#' slotNames(x)
NULL



#' Markdown Links to Results Files
#'
#' Enables looping of results contrast file links for R Markdown.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @return [writeLines()] output.
#' @noRd
.markdownResultsLinks <- function(object, headerLevel = 2L) {
    assert_is_all_of(object, "DESeqResultsTables")
    validObject(object)
    assertIsImplicitInteger(headerLevel)

    markdownHeader(contrast, level = headerLevel, asis = TRUE)

    # Prioritze `dropboxFiles` over `localFiles` for path return.
    # FIXME Need to switch to using slots.
    if ("dropboxFiles" %in% names(object)) {
        # nocov start
        # Using local Dropbox token for code coverage here.
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



.writeResultsTables <- function(
    object,
    dir = ".",
    dropboxDir = NULL,
    rdsToken = rdsToken
) {
    # Write local files to tempdir if Dropbox mode is enabled
    if (is_a_string(dropboxDir)) {
        dir <- tempdir()  # nocov
    } else {
        dir <- initializeDirectory(dir)
    }

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
    # FIXME
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
# TODO By default loop across all results if it's not set here.

#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature("DESeqAnalysis"),
    function(
        object,
        results,
        lfcShrink = TRUE,
        write = FALSE,
        dir = ".",
        dropboxDir = NULL,
        rdsToken = NULL,
        markdown = FALSE,
        headerLevel = 2L,
    ) {
        validObject(object)
        assert_is_a_bool(write)
        assert_is_a_bool(markdown)

        # Extract internal parameters from DESeqResults object -----------------
        results <- .matchResults(
            object = object,
            results = results,
            lfcShrink = lfcShrink
        )
        # Coerce to `DataFrame`.
        # We'll regenerate a modified `DESeqResults` from this below.
        data <- as(results, "DataFrame")

        # Add row annotations --------------------------------------------------
        message("Adding `rowData()` annotations (atomic columns only)")
        rowData <- sanitizeRowData(rowData(object@data))
        # DESeq2 includes additional information in `rowData()` that isn't
        # informative for a user, and doesn't need to be included in the CSV.
        # Use our `bcb_small` example dataset to figure out which columns
        # are worth including.
        keep <- intersect(
            x = colnames(rowData),
            y = colnames(rowData(bcbioRNASeq::bcb_small))
        )
        rowData <- rowData[, keep, drop = FALSE]
        assert_all_are_true(vapply(rowData, is.atomic, logical(1L)))
        assert_is_non_empty(rowData)
        assert_are_disjoint_sets(colnames(data), colnames(rowData))
        data <- cbind(data, rowData)

        # Add the DESeq2 normalized counts -------------------------------------
        message("Adding DESeq2 normalized counts")
        counts <- counts(object@data, normalized = TRUE)
        assert_are_disjoint_sets(colnames(data), colnames(counts))
        assert_are_identical(rownames(data), rownames(counts))
        data <- cbind(data, counts)

        # Regenerate modified DESeqResults -------------------------------------
        # Consider adding information to elementMetadata.
        results <- DESeqResults(
            DataFrame = data,
            priorInfo = priorInfo(results)
        )

        # Return DESeqResultsTables --------------------------------------------
        tables <- resultsTables(results)

        # Write the files to disk in CSV format.
        # Assign back to keep track of file paths (either local or Dropbox).
        if (isTRUE(write)) {
            tables <- .writeResultsTables(
                object = tables,
                dir = dir,
                dropboxDir = dropboxDir,
                rdsToken = rdsToken
            )
        }

        # Include Markdown links, if desired.
        if (isTRUE(markdown)) {
            .markdownResultsLinks(tables, headerLevel = headerLevel)
        }

        tables
    }
)
