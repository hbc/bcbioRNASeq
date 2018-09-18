# FIXME Define SE method for this in basejump.



#' Plot Total Counts
#'
#' @name plotTotalCounts
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotTotalCounts(bcb_small)
NULL



#' @rdname plotTotalCounts
#' @export
setMethod(
    "plotTotalCounts",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups = NULL,
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "total counts"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        counts <- counts(object)
        data <- sampleData(object)
        data[["totalCounts"]] <- colSums(counts)

        p <- ggplot(
                data = as.data.frame(data),
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("totalCounts") / 1e6L,
                    fill = !!sym("interestingGroups")
                )
            ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                y = "counts per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }

        if (isTRUE(flip)) {
            p <- p + coord_flip()
        }

        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(fill = FALSE)
        }

        p
    }
)
