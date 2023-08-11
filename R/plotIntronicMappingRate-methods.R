#' @name plotIntronicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotIntronicMappingRate
#' @note Updated 2022-05-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotIntronicMappingRate(bcb)
NULL



## Updated 2022-05-09.
`plotIntronicMappingRate,bcbioRNASeq` <- # nolint
    function(object,
             interestingGroups = NULL,
             limit = 0.2,
             labels = list(
                 "title" = "Intronic mapping rate",
                 "subtitle" = NULL,
                 "sampleAxis" = NULL,
                 "metricAxis" = "intronic mapping rate (%)"
             ),
             flip = getOption(x = "acid.flip", default = TRUE)) {
        validObject(object)
        assert(
            isProportion(limit),
            isFlag(flip)
        )
        labels <- matchLabels(labels)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = .data[["sampleName"]],
                y = .data[["intronicRate"]] * 100L,
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
            limit <- limit * 100L
            if (limit < 100L) {
                p <- p + acid_geom_abline(yintercept = limit)
            }
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



#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    f = "plotIntronicMappingRate",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotIntronicMappingRate,bcbioRNASeq`
)
