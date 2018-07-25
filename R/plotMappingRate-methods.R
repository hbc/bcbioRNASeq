#' Plot Mapping Rate
#'
#' The genomic mapping rate represents the percentage of reads mapping to the
#' reference genome. Low mapping rates are indicative of sample contamination,
#' poor sequencing quality or other artifacts.
#'
#' @name plotMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotMappingRate(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotMappingRate
#' @export
setMethod(
    "plotMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        limit = 90L,
        fill = NULL,
        flip = TRUE,
        title = "mapping rate"
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / !!sym("totalReads") * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            ylim(0L, 100L) +
            labs(
                title = title,
                x = NULL,
                y = "mapping rate (%)",
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
