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



.topTable <- function(
    object,
    direction = c("up", "down"),
    n = 50L
) {
    assert_is_all_of(object, "DESeqResultsTables")
    validObject(object)
    direction <- match.arg(direction)
    assertIsImplicitInteger(n)

    data <- slot(object, name = camel(paste("deg", direction)))
    assert_is_all_of(data, "DataFrame")

    # Early return `NULL` when there are no significant DEGs
    if (!nrow(data)) {
        return(invisible())  # nocov
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



.topTables.DESeqAnalysis <-  # nolint
    function(
        object,
        results = 1L,
        lfcShrink = TRUE,
        n = 50L
    ) {
        results <- .matchResults(
            object = object,
            results = results,
            lfcShrink = lfcShrink
        )
        geneCols <- c("geneName", "geneBiotype", "description")
        assert_are_disjoint_sets(geneCols, colnames(results))

        # Add useful gene annotations.
        rowData <- rowData(object@data)
        assert_is_subset(geneCols, colnames(rowData))

        rowData <- rowData[, geneCols, drop = FALSE]
        # Don't use cbind or merge...it will coerce to `DataFrame`.
        results[["geneName"]] <- rowData[["geneName"]]
        results[["geneBiotype"]] <- rowData[["geneBiotype"]]
        results[["description"]] <- rowData[["description"]]
        assert_is_all_of(results, "DESeqResults")

        # Using DESeqResults method.
        do.call(
            what = topTables,
            args = list(
                object = DESeqResultsTables(results),
                n = n
            )
        )
    }



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



.topTables.DESeqResultsTables <-  # nolint
    function(object, n = 50L) {
        validObject(object)
        assertIsAnImplicitInteger(n)
        contrast <- contrastName(object)
        assert_is_a_string(contrast)

        # Upregulated.
        up <- .topTable(object, direction = "up", n = n)
        if (length(up)) {
            show(kable(
                x = as.data.frame(up),
                caption = paste(contrast, "(upregulated)")
            ))
        }

        # Downregulated.
        down <- .topTable(object, direction = "down", n = n)
        if (length(down)) {
            show(kable(
                x = as.data.frame(down),
                caption = paste(contrast, "(downregulated)")
            ))
        }

        # Invisibly return list containing the subsets.
        invisible(list(up = up, down = down))
    }



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("DESeqAnalysis"),
    definition = .topTables.DESeqAnalysis
)



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("DESeqResults"),
    definition = .topTables.DESeqResults
)



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqResultsTables"),
    definition = .topTables.DESeqResultsTables
)
