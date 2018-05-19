#' Top Tables of Differential Expression Results
#'
#' @name topTables
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param object [resultsTables()] return `list`.
#' @param n Number genes to report.
#' @param coding Whether to only return coding genes.
#'
#' @return `kable`.
#'
#' @examples
#' # DESeqResults ====
#' # Minimal return
#' topTables(res_small, n = 5L)
#'
#' # resultsTables list ====
#' # Return with gene annotations and DESeq2 normalized counts
#' x <- resultsTables(
#'     results = res_small,
#'     counts = dds_small
#' )
#' topTables(x, n = 5L)
NULL



# Constructors =================================================================
.subsetTop <- function(
    object,
    n = 50L,
    coding = FALSE
) {
    assert_is_data.frame(object)
    assert_has_colnames(object)
    assert_has_rows(object)
    assertIsImplicitInteger(n)
    assert_is_a_bool(coding)
    # Note that `geneName` and `description` columns are optional
    requiredCols <- c("geneID", "baseMean", "log2FoldChange", "padj")
    assert_is_subset(requiredCols, colnames(object))

    if (isTRUE(coding)) {
        assert_is_subset("broadClass", colnames(object))
        object <- object %>%
            .[.[["broadClass"]] == "coding", , drop = FALSE]
    }

    # Early return NULL when there are no significant DEGs
    if (!nrow(object)) {
        return(NULL)
    }

    keepCols <- c(requiredCols, c("geneName", "geneBiotype", "description"))
    data <- object %>%
        as_tibble() %>%
        remove_rownames() %>%
        head(n = n) %>%
        mutate(
            baseMean = round(!!sym("baseMean")),
            log2FoldChange = format(!!sym("log2FoldChange"), digits = 3L),
            padj = format(!!sym("padj"), digits = 3L, scientific = TRUE)
        ) %>%
        .[, which(colnames(.) %in% keepCols)] %>%
        # Shorten `log2FoldChange` to `lfc` to keep column width compact
        rename(lfc = !!sym("log2FoldChange"))

    # Sanitize the description, if necessary
    if ("description" %in% colnames(data)) {
        # Remove symbol information in description, if present
        data[["description"]] <- gsub(
            pattern = " \\[.+\\]$",
            replacement = "",
            x = data[["description"]]
        )
    }

    # Ensure `gene*` columns appear first
    data <- select(data, starts_with("gene"), everything())

    data
}



# Methods ======================================================================
#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqResults"),
    function(
        object,
        n = 50L,
        coding = FALSE
    ) {
        contrast <- contrastName(object)
        padj <- object %>%
            as.data.frame() %>%
            rownames_to_column("geneID") %>%
            # Remove any rows with NA P values
            .[complete.cases(.), ] %>%
            .[order(.[["padj"]]), ]
        up <- padj %>%
            .[.[["log2FoldChange"]] > 0L, , drop = FALSE] %>%
            .subsetTop(n = n, coding = coding)
        down <- padj %>%
            .[.[["log2FoldChange"]] < 0L, , drop = FALSE] %>%
            .subsetTop(n = n, coding = coding)
        if (!is.null(up)) {
            show(kable(
                up,
                caption = paste(contrast, "(upregulated)")
            ))
        }
        if (!is.null(down)) {
            show(kable(
                down,
                caption = paste(contrast, "(downregulated)")
            ))
        }
    }
)



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("list"),
    function(
        object,
        n = 50L,
        coding = FALSE
    ) {
        assert_is_list(object)
        assert_is_subset(
            c("all", "deg", "degLFCDown", "degLFCUp"),
            names(object)
        )
        up <- .subsetTop(
            object[["degLFCUp"]],
            n = n,
            coding = coding
        )
        down <- .subsetTop(
            object[["degLFCDown"]],
            n = n,
            coding = coding
        )
        contrast <- object[["contrast"]]
        if (!is.null(up)) {
            show(kable(
                up,
                caption = paste(contrast, "(upregulated)")
            ))
        }
        if (!is.null(down)) {
            show(kable(
                down,
                caption = paste(contrast, "(downregulated)")
            ))
        }
    }
)
