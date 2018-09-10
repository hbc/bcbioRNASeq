#' Top Tables of Differential Expression Results
#'
#' @name topTables
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param n `scalar integer`. Number genes to report.
#'
#' @return `kable`.
#'
#' @examples
#' # DESeqAnalysis ====
#' topTables(deseq_small, results = 1L, n = 5L)
#'
#' # DESeqResults ====
#' object <- deseq_small@results[[1L]]
#' topTables(object, n = 5L)
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
        "description",
        "geneName",
        "geneBiotype",
        "rowname"
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



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqAnalysis"),
    function(
        object,
        results,
        n = 50L
    ) {
        results <- .matchResults(object, results)
        do.call(
            what = topTables,
            args = list(
                object = resultsTables(results),
                n = n
            )
        )
    }
)



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqResults"),
    function(object, n = 50L) {
        do.call(
            what = topTables,
            args = list(
                object = resultsTables(object),
                n = n
            )
        )
    }
)



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqResultsTables"),
    function(object, n = 50L) {
        validObject(object)
        assertIsAnImplicitInteger(n)
        contrast <- contrastName(object)
        assert_is_a_string(contrast)

        # Upregulated
        up <- .topTable(object, direction = "up", n = n)
        if (length(up)) {
            show(kable(
                x = as.data.frame(up),
                caption = paste(contrast, "(upregulated)")
            ))
        }

        # Downregulated
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
)
