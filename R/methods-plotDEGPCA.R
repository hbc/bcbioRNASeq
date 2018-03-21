# TODO Check for interestingGroups stashed in DESeqTransform

#' Plot DEG PCA
#'
#' @name plotDEGPCA
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inherit plotPCA
#'
#' @examples
#' # DESeqResults, DESeqTransform ====
#' plotDEGPCA(
#'     object = res_small,
#'     counts = rld_small,
#'     interestingGroups = interestingGroups(bcb_small),
#'     label = TRUE
#' )
NULL



# Methods ======================================================================
#' @rdname plotDEGPCA
#' @export
setMethod(
    "plotDEGPCA",
    signature("DESeqResults"),
    function(
        object,
        counts,
        interestingGroups = "sampleName",
        lfc = 0L,
        color = scale_color_viridis(discrete = TRUE),
        label = FALSE,
        title = "deg pca",
        return = c("ggplot", "data.frame")
    ) {
        # Passthrough: interestingGroups, color, label, title
        assert_is_all_of(counts, "DESeqTransform")
        assert_is_a_number(lfc)
        assert_all_are_non_negative(lfc)
        return <- match.arg(return)

        # Get the DE gene vector using `resultsTables()`
        list <- resultsTables(
            object,
            lfc = lfc,
            rowData = NULL,
            summary = FALSE,
            write = FALSE
        )
        genes <- c(
            list[["degLFCUp"]][["geneID"]],
            list[["degLFCDown"]][["geneID"]]
        )

        # Subset the matrix
        counts <- counts[genes, , drop = FALSE]

        data <- plotPCA(
            object = counts,
            intgroup = interestingGroups,
            ntop = nrow(counts),
            returnData = TRUE
        ) %>%
            camel()

        if (return == "data.frame") {
            return(data)
        }

        # Use `sampleName` for plot labels
        if (isTRUE(label)) {
            assert_is_subset("sampleName", colnames(colData(counts)))
            data[["label"]] <- colData(counts)[, "sampleName", drop = TRUE]
        }

        .plotPCA.ggplot(
            object = data,
            color = color,
            label = label,
            group = interestingGroups,
            title = title
        )
    }
)
