#' @name topTables
#' @inherit basejump::topTables
#' @author Michael Steinbaugh
#'
#' @inheritParams params
#' @param n `scalar integer`. Number of genes (per direction) to report.
#'
#' @examples
#' data(deseq)
#' topTables(deseq, results = 1L, n = 5L)
NULL



#' @importFrom basejump topTables
#' @aliases NULL
#' @export
basejump::topTables



# DESeqResults =================================================================
topTables.DESeqResults <-  # nolint
    function(object, n = 10L) {
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
    definition = topTables.DESeqResults
)



# DESeqResultsTables ===========================================================
.topTable <-  # nolint
    function(
        object,
        direction = c("up", "down"),
        n
    ) {
        assert_that(is(object, "DESeqResultsTables"))
        validObject(object)
        direction <- match.arg(direction)
        assertIsImplicitInteger(n)

        results <- slot(object, "results")
        assert_is_all_of(results, "DESeqResults")

        # Get the differentially expressed genes.
        deg <- slot(object, "deg")[[direction]]
        # Early return when there are no significant DEGs.
        if (!has_length(deg)) {
            warning(paste0(
                "No significant ", direction, "-regulated genes."
            ), call. = FALSE)
            return(invisible())
        }

        # Coerce DESeqResults to DEG subset DataFrame.
        data <- as(results, "DataFrame")
        keep <- c("baseMean", "log2FoldChange", "padj")
        assert_is_subset(keep, colnames(data))
        assert_is_subset(deg, rownames(data))
        data <- data[deg, keep, drop = FALSE]

        # Join helpful row annotations, if defined.
        rowRanges <- slot(object, "rowRanges")
        # Note that `rlang::has_length()` doesn't return FALSE here.
        if (length(rowRanges) > 0L) {
            rowData <- as(rowRanges, "DataFrame")
            # Consider using broadClass here instead of geneBiotype.
            keep <- intersect(
                x = c("geneName", "geneBiotype", "description"
                ),
                y = colnames(rowData)
            )
            assert_is_subset(rownames(data), rownames(rowData))
            rowData <- rowData[rownames(data), keep, drop = FALSE]
            if ("description" %in% colnames(rowData)) {
                # Remove symbol information in brackets.
                rowData[["description"]] <- gsub(
                    pattern = " \\[.+\\]$",
                    replacement = "",
                    x = rowData[["description"]]
                )
                # Truncate to max 50 characters.
                rowData[["description"]] <- str_trunc(
                    string = as.character(rowData[["description"]]),
                    width = 50L,
                    side = "right"
                )
            }
            data <- cbind(data, rowData)
        }

        # Modify using dplyr and return as DataFrame.
        data %>%
            as("tbl_df") %>%
            head(n = n) %>%
            mutate(
                baseMean = round(!!sym("baseMean")),
                log2FoldChange = format(!!sym("log2FoldChange"), digits = 3L),
                padj = format(!!sym("padj"), digits = 3L, scientific = TRUE)
            ) %>%
            # Shorten `log2FoldChange` to `lfc` to keep column width compact.
            rename(lfc = !!sym("log2FoldChange")) %>%
            # Ensure `gene*` annotation columns appear first.
            select(starts_with("gene"), everything()) %>%
            # This coercion will automatically set rownames.
            as("DataFrame")
    }
formals(.topTable)[["n"]] <-
    formals(topTables.DESeqResults)[["n"]]



topTables.DESeqResultsTables <-  # nolint
    function(object, n) {
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
formals(topTables.DESeqResultsTables)[["n"]] <-
    formals(topTables.DESeqResults)[["n"]]



#' @rdname topTables
#' @export
setMethod(
    "topTables",
    signature("DESeqResultsTables"),
    definition = topTables.DESeqResultsTables
)



# DESeqAnalysis ================================================================
topTables.DESeqAnalysis <-  # nolint
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
                    lfcShrink = lfcShrink
                ),
                n = n
            )
        )
    }
formals(topTables.DESeqAnalysis)[["n"]] <-
    formals(topTables.DESeqResults)[["n"]]



#' @rdname topTables
#' @export
setMethod(
    f = "topTables",
    signature = signature("DESeqAnalysis"),
    definition = topTables.DESeqAnalysis
)
