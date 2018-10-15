# FIXME Need to fix these for DESeqResultsTables update.



#' Top Tables of Differential Expression Results
#'
#' @name topTables
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @param n `scalar integer`. Number of genes (per direction) to report.
#'
#' @return `kable`. Markdown tables.
#'
#' @examples
#' data(deseq_small)
#' topTables(deseq_small, results = 1L, n = 5L)
NULL



# DESeqResults =================================================================
.topTables.DESeqResults <-  # nolint
    function(object, n = 50L) {
        do.call(
            what = topTables,
            args = list(
                object = DESeqResultsTables(object),
                n = n
            )
        )
    }



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("DESeqResults"),
    definition = .topTables.DESeqResults
)



# DESeqResultsTables ===========================================================
# FIXME Improve support for joining rowData.

# geneCols <- c("geneName", "geneBiotype", "description")
# assert_are_disjoint_sets(geneCols, colnames(results))
#
# # Add useful gene annotations.
# rowData <- rowData(object@data)
# assert_is_subset(geneCols, colnames(rowData))
#
# rowData <- rowData[, geneCols, drop = FALSE]
# # Don't use cbind or merge...it will coerce to `DataFrame`.
# results[["geneName"]] <- rowData[["geneName"]]
# results[["geneBiotype"]] <- rowData[["geneBiotype"]]
# results[["description"]] <- rowData[["description"]]
# assert_is_all_of(results, "DESeqResults")

.topTable.DESeqResultsTables <-  # nolint
    function(
        object,
        direction = c("up", "down"),
        n
    ) {
        stopifnot(is(object, "DESeqResultsTables"))
        validObject(object)
        direction <- match.arg(direction)
        assertIsImplicitInteger(n)

        results <- object@results
        assert_is_all_of(results, "DESeqResults")

        rownames <- object@deg[[direction]]
        assert_is_subset(rownames, rownames(results))
        data <- as(results, "DataFrame")[rownames, , drop = FALSE]

        # Early return when there are no significant DEGs.
        if (!nrow(data)) {
            warning(paste0(
                "No significant ", direction, "-regulated genes."
            ), call. = FALSE)
            return(invisible())
        }

        requiredCols <- c(
            "baseMean",
            "log2FoldChange",
            "padj"
        )
        keepCols <- c(
            requiredCols,
            "rowname",
            "geneName",
            "geneBiotype",
            "description"
        )
        assert_is_subset(requiredCols, colnames(data))

        # Sanitize the description, if necessary.
        if ("description" %in% colnames(data)) {
            # Remove symbol information in description, if present
            data[["description"]] <- gsub(
                pattern = " \\[.+\\]$",
                replacement = "",
                x = data[["description"]]
            )
        }

        # Coerce to `tbl_df` and use dplyr functions, then return `data.frame`.
        data %>%
            as("tbl_df") %>%
            head(n = n) %>%
            mutate(
                baseMean = round(!!sym("baseMean")),
                log2FoldChange = format(!!sym("log2FoldChange"), digits = 3L),
                padj = format(!!sym("padj"), digits = 3L, scientific = TRUE)
            ) %>%
            # Select only the desired keep columns (see above).
            .[, which(colnames(.) %in% keepCols)] %>%
            # Shorten `log2FoldChange` to `lfc` to keep column width compact.
            rename(lfc = !!sym("log2FoldChange")) %>%
            # Ensure `gene*` annotation columns appear first.
            select(starts_with("gene"), everything()) %>%
            # This coercion will automatically set rownames.
            as("DataFrame")
    }
formals(.topTable.DESeqResultsTables)[["n"]] <-
    formals(.topTables.DESeqResults)[["n"]]



.topTables.DESeqResultsTables <-  # nolint
    function(object, n) {
        validObject(object)
        assertIsAnImplicitInteger(n)
        contrast <- contrastName(object)
        assert_is_a_string(contrast)

        # Upregulated.
        up <- .topTable.DESeqResultsTables(
            object = object,
            direction = "up",
            n = n
        )
        if (length(up)) {
            show(kable(
                x = as.data.frame(up),
                caption = paste(contrast, "(upregulated)")
            ))
        }

        # Downregulated.
        down <- .topTable.DESeqResultsTables(
            object = object,
            direction = "down",
            n = n
        )
        if (length(down)) {
            show(kable(
                x = as.data.frame(down),
                caption = paste(contrast, "(downregulated)")
            ))
        }

        # Invisibly return list containing the subsets.
        invisible(list(up = up, down = down))
    }
formals(.topTables.DESeqResultsTables)[["n"]] <-
    formals(.topTables.DESeqResults)[["n"]]



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqResultsTables"),
    definition = .topTables.DESeqResultsTables
)



# DESeqAnalysis ================================================================
.topTables.DESeqAnalysis <-  # nolint
    function(
        object,
        results = 1L,
        lfcShrink = TRUE,
        n
    ) {
        # Using DESeqResultsTables method.
        do.call(
            what = topTables,
            args = list(
                object = DESeqResultsTables(
                    object = object,
                    results = results,
                    lfcShrink = lfcShrink,
                    rowData = TRUE,
                    counts = FALSE
                ),
                n = n
            )
        )
    }
formals(.topTables.DESeqAnalysis)[["n"]] <-
    formals(.topTables.DESeqResults)[["n"]]



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("DESeqAnalysis"),
    definition = .topTables.DESeqAnalysis
)
