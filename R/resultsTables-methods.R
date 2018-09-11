# TODO Consider making `resultsTables()` have fewer options...split these out
# into separate functions.



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
#' @param rowData `boolean`. Include gene annotations.
#' @param counts `boolean`. Include DESeq2 normalized counts.
#' @param write `boolean`. Write CSV files to disk.
#' @param dropboxDir `string` or `NULL`. Dropbox directory path where to archive
#'   the results tables for permanent storage (e.g. Stem Cell Commons). When
#'   this option is enabled, unique links per file are generated internally with
#'   the rdrop2 package. Note that local files are written to [base::tempdir()]
#'   and the `dir` argument is ignored, if this is enabled.
#' @param rdsToken `string` or `NULL`. RDS file token to use for Dropbox
#'   authentication. If set `NULL` and `dropboxDir` is defined, then an
#'   interactive prompt will appear requesting authorization.
#' @param markdown `boolean`. Include Markdown links to file paths. Only applies
#'   when `write = TRUE`.
#'
#' @return `DESeqResultsTables`.
#'
#' @examples
#' # DESeqAnalysis ====
#' # This is the recommended default method.
#' x <- resultsTables(res)
#'
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
.markdownResultsTables <- function(
    object,
    headerLevel = 2L
) {
    assert_is_all_of(object, "DESeqResultsTables")
    validObject(object)
    assertIsImplicitInteger(headerLevel)

    # Include a contrast header, which is useful for looping.
    contrast <- contrastName(object)
    markdownHeader(contrast, level = headerLevel, asis = TRUE)
    headerLevel <- headerLevel + 1L

    # CSV file links -----------------------------------------------------------
    # Prioritze `dropboxFiles` over `localFiles`.
    if (length(object@dropboxFiles) > 0L) {
        # nocov start
        # Using local Dropbox token for code coverage here.
        paths <- vapply(
            X = object@dropboxFiles,
            FUN = function(x) {
                x[["url"]]
            },
            FUN.VALUE = "character"
        )
        basenames <- gsub("\\?.*$", "", basename(paths))
        # nocov end
    } else if (length(object@localFiles) > 0L) {
        paths <- object@localFiles
        basenames <- basename(paths)
    } else {
        stop("Object doesn't contain saved file paths")
    }
    names(basenames) <- names(paths)

    markdownHeader("CSV file links", level = headerLevel, asis = TRUE)
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
            "(", paths[["degUp"]], "): ",
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
    rdsToken = NULL
) {
    validObject(object)
    assertIsAStringOrNULL(dropboxDir)
    # Write local files to tempdir if Dropbox mode is enabled.
    if (is_a_string(dropboxDir)) {
        dir <- tempdir()  # nocov
    } else {
        dir <- initializeDirectory(dir)
    }

    # Extract the results tables from the object.
    tables <- flatFiles(object)[c("all", "deg", "degUp", "degDown")]

    # Local files (required) -------------------------------------------
    stem <- snake(contrastName(object))
    localFiles <- file.path(
        dir,
        paste0(stem, "_", snake(names(tables)), ".csv.gz")
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

    invisible(mapply(
        x = tables,
        path = localFiles,
        FUN = function(x, path) {
            x <- as(x, "tbl_df")
            assert_is_subset("rowname", colnames(x))
            write_csv(x = x, path = path)
        },
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    ))

    # Check that all writes were successful.
    assert_all_are_existing_files(localFiles)

    # Assign the local file paths to the object.
    # Slot these only when we're not saving to Dropbox.
    if (is.null(dropboxDir)) {
        slot(object, "localFiles") <- localFiles
    }

    # Copy to Dropbox (optional) ---------------------------------------
    if (is.character(dropboxDir)) {
        # nocov start
        # Using a local token to cover this code.
        dropboxFiles <- copyToDropbox(
            files = localFiles,
            dir = dropboxDir,
            rdsToken = rdsToken
        )
        assert_is_list(dropboxFiles)
        slot(object, "dropboxFiles") <- dropboxFiles
        # nocov end
    }

    object
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



#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature("DESeqAnalysis"),
    function(
        object,
        results,
        lfcShrink = TRUE,
        rowData = TRUE,
        counts = TRUE,
        write = FALSE,
        dir = ".",
        dropboxDir = NULL,
        rdsToken = NULL,
        markdown = FALSE,
        headerLevel = 2L
    ) {
        validObject(object)
        assert_is_a_bool(rowData)
        assert_is_a_bool(counts)
        assert_is_a_bool(write)
        assert_is_a_bool(markdown)

        # Match the DESeqResults object.
        results <- .matchResults(
            object = object,
            results = results,
            lfcShrink = lfcShrink
        )

        # Add columns (optional) -----------------------------------------------
        if (
            isTRUE(rowData) ||
            isTRUE(counts)
        ) {
            # Coerce to `DataFrame`.
            # We'll regenerate a modified `DESeqResults` from this below.
            data <- as(results, "DataFrame")

            # Row annotations
            if (isTRUE(rowData)) {
                message("Adding `rowData()` annotations (atomic columns only)")
                rowData <- sanitizeRowData(rowData(object@data))
                # DESeq2 includes additional information in `rowData()` that
                # isn't informative for a user, and doesn't need to be included
                # in the CSV. Use our `bcb_small` example dataset to figure out
                # which columns are worth including.
                keep <- intersect(
                    x = colnames(rowData),
                    y = colnames(rowData(bcbioRNASeq::bcb_small))
                )
                rowData <- rowData[, keep, drop = FALSE]
                assert_all_are_true(vapply(rowData, is.atomic, logical(1L)))
                assert_is_non_empty(rowData)
                assert_are_disjoint_sets(colnames(data), colnames(rowData))
                data <- cbind(data, rowData)
            }

            # DESeq2 normalized counts
            if (isTRUE(counts)) {
                message("Adding DESeq2 normalized counts")
                counts <- counts(object@data, normalized = TRUE)
                assert_are_disjoint_sets(colnames(data), colnames(counts))
                assert_are_identical(rownames(data), rownames(counts))
                data <- cbind(data, counts)
            }

            # Regenerate the DESeqResults.
            # Consider adding information to elementMetadata?
            results <- DESeqResults(
                DataFrame = data,
                priorInfo = priorInfo(results)
            )
        }

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
            # Include Markdown links, if desired.
            # FIXME Break this out into a separate function.
            if (isTRUE(markdown)) {
                .markdownResultsTables(
                    object = tables,
                    headerLevel = headerLevel
                )
            }
        }

        tables
    }
)
