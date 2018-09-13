#' Plot Counts Per Gene
#'
#' Generally, we expect similar count spreads for all genes between samples
#' unless the library sizes or total RNA expression are different.
#'
#' @section TMM:
#' We recommend visualizing counts normalized with the Trimmed Mean of M-Values
#' (TMM) method here (Robinson, et al., 2010). TMM normalization equates the
#' overall expression levels of genes between samples under the assumption that
#' the majority of them are not differentially expressed. Therefore, by
#' normalizing for total RNA expression by sample, we expect the spread of the
#' TMM-normalized counts per gene to be similar for every sample.
#'
#' @name plotCountsPerGene
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @export
#'
#' @inheritParams general
#' @param geom `string`. Type of ggplot2 geometric object to use.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotCountsPerGene(bcb_small)
NULL



#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups = NULL,
        normalized = c("tmm", "vst", "rlog", "tpm", "rle"),
        geom = c("violin", "density", "boxplot"),
        color = getOption("bcbio.discrete.color", NULL),
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "counts per gene"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        normalized <- match.arg(normalized)
        geom <- match.arg(geom)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        data <- .meltCounts(object = object, normalized = normalized)
        count <- length(unique(data[["rowname"]]))

        # Subtitle
        if (is_a_string(title)) {
            subtitle <- paste(count, "non-zero genes")
        } else {
            subtitle <- NULL
        }

        # Construct the ggplot.
        p <- ggplot(data = data)
        countsAxisLabel <- paste(normalized, "counts (log2)")

        if (geom == "violin") {
            p <- p +
                geom_violin(
                    mapping = aes(
                        x = !!sym("sampleName"),
                        y = !!sym("counts"),
                        fill = !!sym("interestingGroups")
                    ),
                    color = "black",
                    scale = "width"
                ) +
                labs(x = NULL, y = countsAxisLabel)
        } else if (geom == "density") {
            p <- p +
                geom_density(
                    mapping = aes(
                        x = !!sym("counts"),
                        group = !!sym("interestingGroups"),
                        color = !!sym("interestingGroups")
                    ),
                    fill = NA,
                    size = 1L
                ) +
                labs(x = countsAxisLabel)
        } else if (geom == "boxplot") {
            p <- p +
                geom_boxplot(
                    mapping = aes(
                        x = !!sym("sampleName"),
                        y = !!sym("counts"),
                        fill = !!sym("interestingGroups")
                    ),
                    color = "black",
                    outlier.shape = NA
                ) +
                labs(x = countsAxisLabel, y = NULL)
        }

        # Add the axis and legend labels.
        p <- p +
            labs(
                title = title,
                subtitle = subtitle,
                color = paste(interestingGroups, collapse = ":\n"),
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }

        # Flip the axis for plots with counts on y-axis, if desired.
        if (isTRUE(flip) && !geom %in% "density") {
            p <- p + coord_flip()
        }

        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(color = FALSE, fill = FALSE)
        }

        p
    }
)
