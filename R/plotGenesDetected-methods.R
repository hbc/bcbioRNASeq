#' Plot Genes Detected
#'
#' @name plotGenesDetected
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotGenesDetected(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        limit = 0L,
        minCounts = 1L,
        fill = getOption("bcbio.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "genes detected"
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        assert_is_character(interestingGroups)
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsAnImplicitInteger(minCounts)
        assert_all_are_in_range(minCounts, lower = 1L, upper = Inf)
        assert_all_are_non_negative(minCounts)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        counts <- counts(object, normalized = FALSE)
        geneCount <- colSums(counts >= minCounts)

        p <- metrics(object) %>%
            mutate(geneCount = !!geneCount) %>%
            ggplot(
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("geneCount"),
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
                y = "gene count",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            p <- p + bcbio_geom_abline(yintercept = limit)
        }

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
