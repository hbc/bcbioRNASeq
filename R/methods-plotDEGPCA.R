#' Plot DEG PCA
#'
#' @name plotDEGPCA
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotPCA
#' @inheritParams plotDEGHeatmap
#'
#' @examples
#' # DESeqResults, DESeqTransform ====
#' plotDEGPCA(
#'     object = res_small,
#'     counts = rld_small,
#'     interestingGroups = interestingGroups(bcb_small)
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
        return = c("ggplot", "data.frame")
    ) {
        # Passthrough: interestingGroups, color, label, return
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

        if (return == "data.frame") {
            returnData <- TRUE
        } else {
            returnData <- FALSE
        }

        # TODO Update to use our PCA code
        plotPCA(
            object = counts,
            intgroup = interestingGroups,
            ntop = nrow(counts),
            returnData = returnData
        )
    }
)
