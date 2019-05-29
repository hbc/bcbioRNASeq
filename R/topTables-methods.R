#' @name topTables
#' @inherit bioverbs::topTables
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @param object `DESeqResults` or `list`. For `list` method, must use
#'   [resultsTables()] return.
#' @param n `scalar integer`. Number genes to report.
#' @param coding `boolean`. Whether to only return coding genes.
#'
#' @return `kable`.
#'
#' @examples
#' # DESeqResults ====
#' # Minimal return
#' topTables(res_small, n = 5L)
#'
#' # `resultsTables()` return list ====
#' # Return with gene annotations and DESeq2 normalized counts
#' res_tbl <- resultsTables(
#'     results = res_small,
#'     counts = dds_small
#' )
#' topTables(res_tbl, n = 5L)
NULL



#' @rdname topTables
#' @name topTables
#' @importFrom bioverbs topTables
#' @usage topTables(object, ...)
#' @export
NULL



.subsetTop <- function(
    object,
    direction = c("up", "down"),
    n = 50L,
    coding = FALSE
) {
    data <- as.data.frame(object)
    assert_has_colnames(data)
    assert_has_rows(data)
    assertIsImplicitInteger(n)
    assert_is_a_bool(coding)
    # Note that `geneName` and `description` columns are optional
    requiredCols <- c("geneID", "baseMean", "log2FoldChange", "padj")
    if (!"geneID" %in% colnames(data)) {
        data <- rownames_to_column(data, "geneID")
    }
    assert_is_subset(requiredCols, colnames(data))

    # Coerce to tibble and arrange by adjusted P value
    data <- data %>%
        as_tibble() %>%
        remove_rownames() %>%
        # Remove rows containing NA P value
        filter(!is.na(!!sym("padj"))) %>%
        arrange(!!sym("padj"))

    # Apply direction
    if (direction == "up") {
        data <- filter(data, !!sym("log2FoldChange") > 0L)
    } else if (direction == "down") {
        data <- filter(data, !!sym("log2FoldChange") < 0L)
    }

    # Coding genes only, if desired
    if (isTRUE(coding)) {
        assert_is_subset("broadClass", colnames(data))
        data <- filter(data, !!sym("broadClass") == "coding")
    }

    # Early return NULL when there are no significant DEGs
    if (!nrow(data)) {
        return(NULL)  # nocov
    }

    keepCols <- c(requiredCols, c("geneName", "geneBiotype", "description"))
    data <- data %>%
        head(n = n) %>%
        mutate(
            baseMean = round(!!sym("baseMean")),
            log2FoldChange = format(!!sym("log2FoldChange"), digits = 3L),
            padj = format(!!sym("padj"), digits = 3L, scientific = TRUE)
        ) %>%
        .[, which(colnames(.) %in% keepCols)] %>%
        # Shorten `log2FoldChange` to `lfc` to keep column width compact
        rename(lfc = !!sym("log2FoldChange")) %>%
        # Ensure `gene*` columns appear first
        select(starts_with("gene"), everything())

    # Sanitize the description, if necessary
    if ("description" %in% colnames(data)) {
        # Remove symbol information in description, if present
        data[["description"]] <- gsub(
            pattern = " \\[.+\\]$",
            replacement = "",
            x = data[["description"]]
        )
    }

    data
}



topTables.DESeqResults <-  # nolint
    function(
        object,
        n = 50L,
        coding = FALSE
    ) {
        assertIsAnImplicitInteger(n)
        assert_is_a_bool(coding)
        contrast <- contrastName(object)

        up <- .subsetTop(object, direction = "up", n = n, coding = coding)
        down <- .subsetTop(object, direction = "down", n = n, coding = coding)

        if (length(up)) {
            show(kable(
                x = up,
                caption = paste(contrast, "(upregulated)")
            ))
        }
        if (length(down)) {
            show(kable(
                x = down,
                caption = paste(contrast, "(downregulated)")
            ))
        }

        # Invisibly return list containing the subset data.frames
        invisible(list(up = up, down = down))
    }



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("DESeqResults"),
    definition = topTables.DESeqResults
)



topTables.list <-  # nolint
    function(
        object,
        n = 50L,
        coding = FALSE
    ) {
        assertIsAnImplicitInteger(n)
        assert_is_a_bool(coding)
        assert_is_subset(
            x = c("all", "deg", "degLFCDown", "degLFCUp"),
            y = names(object)
        )
        contrast <- object[["contrast"]]
        assert_is_a_string(contrast)

        up <- .subsetTop(
            object = object[["degLFCUp"]],
            direction = "up",
            n = n,
            coding = coding
        )
        down <- .subsetTop(
            object = object[["degLFCDown"]],
            direction = "down",
            n = n,
            coding = coding
        )

        if (length(up)) {
            show(kable(
                x = up,
                caption = paste(contrast, "(upregulated)")
            ))
        }
        if (length(down)) {
            show(kable(
                x = down,
                caption = paste(contrast, "(downregulated)")
            ))
        }

        # Invisibly return list containing the subset data.frames
        invisible(list(up = up, down = down))
    }



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("list"),
    definition = topTables.list
)
