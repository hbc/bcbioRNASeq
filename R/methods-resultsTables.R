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

    all <- results %>%
        as.data.frame() %>%
        camel()

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
        "deg" = deg,
        "degLFC" = degLFC,
        "degLFCUp" = degLFCUp,
        "degLFCDown" = degLFCDown,
        "all" = all,
        "contrast" = contrast,
        "alpha" = alpha,
        "lfcThreshold" = lfcThreshold
    )
}



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
# Minimal method that is used by other functions inside the package
#' @rdname resultsTables
#' @export
setMethod(
    "resultsTables",
    signature(
        results = "DESeqResults",
        counts = "missingOrNULL"
    ),
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
)



# Only dispatch with advanced params if DESeqDataSet is defined
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
        dir <- initializeDirectory(dir)

        # Extract internal parameters from DESeqResults object =================
        if (missing(alpha)) {
            alpha <- metadata(results)[["alpha"]]
        }
        assert_is_a_number(alpha)
        contrast <- contrastName(results)
        fileStem <- snake(contrast)

        # Now safe to coerce to DataFrame
        results <- as.data.frame(results) %>%
            rownames_to_column("geneID")

        # Add gene annotations (rowData) =======================================
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
        # Coerce to `data.frame` here first instead of `DataFrame` so the
        # GRanges gets collapsed into columns. `DataFrame` nests this data in a
        # column named `X` otherwise, and that won't write to disk in CSV
        # format.
        rowData <- rowRanges %>%
            as.data.frame() %>%
            as("DataFrame")
        # Drop columns already present in results
        rowData <- rowData[, setdiff(colnames(rowData), colnames(results))]

        # Bind annotation columns ==============================================
        results <- cbind(results, rowData)
        results <- cbind(results, counts(counts, normalized = TRUE))

        # Check for overall gene expression with base mean
        baseMeanGt0 <- results %>%
            .[.[["baseMean"]] > 0L, ] %>%
            nrow()
        baseMeanGt1 <- results %>%
            .[.[["baseMean"]] > 1L, ] %>%
            nrow()

        list <- .degList(
            results = results,
            alpha = alpha,
            lfcThreshold = lfcThreshold,
            contrast = contrast
        )

        if (isTRUE(summary)) {
            markdownHeader(
                "Summary statistics",
                level = headerLevel,
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

            # Local files (required) ===========================================
            localFiles <- file.path(
                dir,
                paste0(fileStem, "_", snake(names(tables)), ".csv.gz")
            )
            names(localFiles) <- names(tables)

            # Write the results tables to local directory
            message(paste(
                "Writing", toString(basename(localFiles)), "to", dir
            ))
            mapply(
                x = tables,
                path = localFiles,
                FUN = function(x, path) {
                    assert_are_identical("geneID", colnames(x)[[1L]])
                    write_csv(x = x, path = path)
                }
            )

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
