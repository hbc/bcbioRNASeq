#' Sample PCA Plot for Transformed Data
#'
#' Wrapper for [DESeq2::plotPCA()] that improves principal component analysis
#' (PCA) sample coloring and labeling.
#'
#' @name plotPCA
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotPCA
#'
#' @inheritParams general
#' @param label Superimpose sample text labels on the plot.
#'
#' @seealso
#' - [DESeq2::plotPCA()].
#' - `getMethod("plotPCA", "DESeqTransform")`
#'
#' @return `ggplot` or `data.frame`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotPCA(bcb_small, label = FALSE)
#' plotPCA(bcb_small, label = TRUE)
#'
#' # Select samples
#' plotPCA(
#'     object = bcb_small,
#'     samples = head(colnames(bcb_small), 4L),
#'     label = TRUE
#' )
NULL



# Constructors =================================================================
.plotPCA.ggplot <- function(  # nolint
    object,
    color = scale_color_viridis(discrete = TRUE),
    label = FALSE,
    group = NULL,
    title = NULL
) {
    assert_is_data.frame(object)
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_bool(label)
    assertIsAStringOrNULL(title)

    percentVar <- round(100L * attr(object, "percentVar"))

    p <- ggplot(
        data = object,
        mapping = aes_string(
            x = "pc1",
            y = "pc2",
            color = "group"
        )
    ) +
        geom_point(size = 4L) +
        coord_fixed() +
        labs(
            title = title,
            x = paste0("pc1: ", percentVar[[1L]], "% variance"),
            y = paste0("pc2: ", percentVar[[2L]], "% variance"),
            color = paste(group, collapse = ":\n")
        )

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    if (isTRUE(label)) {
        p <- p +
            geom_text_repel(
                # Color the labels to match the points
                aes_string(label = "label"),
                color = "black",
                fontface = "bold",
                # Add extra padding around each text label
                box.padding = unit(0.5, "lines"),
                # Add extra padding around each data point
                point.padding = unit(1.5, "lines"),
                # Color of the line segments
                segment.color = "gray",
                # Width of the line segments
                segment.size = 0.5,
                # Draw an arrow from the label to the data point
                arrow = arrow(length = unit(0.01, "npc")),
                # Strength of the repulsion force
                force = 1L,
                show.legend = FALSE
            )
    }

    p
}



# Methods ======================================================================
#' @rdname plotPCA
#' @export
setMethod(
    "plotPCA",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("rlog", "vst", "tmm", "tpm"),
        genes = NULL,
        samples = NULL,
        interestingGroups,
        color = scale_color_viridis(discrete = TRUE),
        label = FALSE,
        title = "pca",
        return = c("ggplot", "data.frame"),
        ...
    ) {
        # Legacy arguments =========================================================
        call <- match.call()
        # censorSamples
        if ("censorSamples" %in% names(call)) {
            warn("`censorSamples` is deprecated in favor of `samples`")
            censorSamples <- call[["censorSamples"]]
            assert_is_subset(censorSamples, colnames(object))
            samples <- setdiff(colnames(object), censorSamples)
        }
        # returnData
        if ("returnData" %in% names(call)) {
            warn("`returnData` is deprecated in favor of `return`")
            returnData <- call[["returnData"]]
            if (isTRUE(returnData)) {
                return <- "data.frame"
            }
        }

        # Assert checks ========================================================
        validObject(object)
        normalized <- match.arg(normalized)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertFormalInterestingGroups(colData(object), interestingGroups)
        assertIsCharacterOrNULL(genes)
        assertIsCharacterOrNULL(samples)
        if (!exists(return, inherits = FALSE)) {
            # exists check is needed for legacy `returnData` match above
            return <- match.arg(return)
        }


        # Prepare counts matrix ================================================
        counts <- counts(object, normalized = normalized)
        colData <- colData(object)

        # Subset samples, if desired
        if (length(samples)) {
            assert_is_subset(samples, colnames(object))
            counts <- counts[, samples, drop = FALSE]
            colData <- colData[samples, , drop = FALSE]
        }

        # Subset genes, if desired
        if (length(genes)) {
            assert_is_subset(genes, rownames(counts))
            counts <- counts[genes, , drop = FALSE]
            # Set ntop to the number of genes requested
            ntop <- length(genes)
            inform(paste(
                "Plotting PCA using", ntop, "genes"
            ))
        } else {
            # Recommended DESeq default of most variable genes
            ntop <- 500L
            inform(paste(
                "Plotting PCA using top", ntop, "most variable genes"
            ))
        }

        # Get PCA data using DESeqTransform method =============================
        se <- SummarizedExperiment(assays = counts, colData = colData)
        dt <- DESeqTransform(se)

        # `group` column is generated by `DESeq2::plotPCA()` for
        # `interestingGroups`
        data <- plotPCA(
            object = dt,
            intgroup = interestingGroups,
            ntop = ntop,
            returnData = TRUE
        ) %>%
            camel()

        if (return == "data.frame") {
            return(data)
        }

        # Use `sampleName` for plot labels
        if (isTRUE(label)) {
            data[["label"]] <- colData[, "sampleName", drop = TRUE]
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
