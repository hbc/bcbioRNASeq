#' @name plotMappedReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotMappedReads
#' @note Updated 2022-05-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotMappedReads(bcb)
NULL



## Updated 2022-05-09.
`plotMappedReads,bcbioRNASeq` <- # nolint
    function(object,
             interestingGroups = NULL,
             limit = 10e6L,
             perMillion = TRUE,
             labels = list(
                 "title" = "Mapped reads",
                 "subtitle" = NULL,
                 "sampleAxis" = NULL,
                 "metricAxis" = "reads"
             ),
             flip = getOption(x = "acid.flip", default = TRUE)) {
        validObject(object)
        assert(
            isInt(limit),
            isNonNegative(limit),
            isFlag(perMillion),
            isFlag(flip)
        )
        labels <- matchLabels(labels)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        ## Convert to per million, if desired.
        if (isTRUE(perMillion)) {
            data[["mappedReads"]] <- data[["mappedReads"]] / 1e6L
            if (isString(labels[["metricAxis"]])) {
                labels[["metricAxis"]] <-
                    paste(labels[["metricAxis"]], "(per million)")
            }
        }
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = .data[["sampleName"]],
                y = .data[["mappedReads"]],
                fill = .data[["interestingGroups"]]
            )
        ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad()
        ## Labels.
        labels[["fill"]] <- paste(interestingGroups, collapse = ":\n")
        names(labels)[names(labels) == "sampleAxis"] <- "x"
        names(labels)[names(labels) == "metricAxis"] <- "y"
        p <- p + do.call(what = labs, args = labels)
        ## Limit.
        if (isPositive(limit)) {
            if (isTRUE(perMillion)) {
                assert(
                    isTRUE(limit >= 1e6L),
                    msg = sprintf(
                        "'%s': %s",
                        limit,
                        "Use absolute value (1e7), not per million (1)."
                    )
                )
                limit <- limit / 1e6L
            }
            p <- p + acid_geom_abline(yintercept = limit)
        }
        ## Color palette.
        p <- p + acid_scale_fill_discrete()
        ## Flip.
        if (isTRUE(flip)) {
            p <- p + acid_discrete_coord_flip()
        }
        ## Hide sample name legend.
        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(fill = "none")
        }
        ## Return.
        p
    }



#' @rdname plotMappedReads
#' @export
setMethod(
    f = "plotMappedReads",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotMappedReads,bcbioRNASeq`
)
