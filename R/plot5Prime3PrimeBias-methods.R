#' @name plot5Prime3PrimeBias
#' @author Michael Steinbaugh
#' @inherit AcidGenerics::plot5Prime3PrimeBias
#' @note Updated 2022-05-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plot5Prime3PrimeBias(bcb)
NULL



## Updated 2022-05-09.
`plot5Prime3PrimeBias,bcbioRNASeq` <- # nolint
    function(object,
             interestingGroups = NULL,
             labels = list(
                 "title" = "5'->3' bias",
                 "subtitle" = NULL,
                 "sampleAxis" = NULL,
                 "metricAxis" = "5'->3' bias"
             ),
             flip = getOption(x = "acid.flip", default = TRUE)) {
        validObject(object)
        assert(isFlag(flip))
        labels <- matchLabels(labels)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = .data[["sampleName"]],
                y = .data[["x5x3Bias"]],
                color = .data[["interestingGroups"]]
            )
        ) +
            geom_point(size = 3L) +
            acid_geom_abline(yintercept = 1L)
        ## Labels.
        labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
        names(labels)[names(labels) == "sampleAxis"] <- "x"
        names(labels)[names(labels) == "metricAxis"] <- "y"
        p <- p + do.call(what = labs, args = labels)
        ## Color palette.
        p <- p + acid_scale_color_discrete()
        ## Flip.
        if (isTRUE(flip)) {
            p <- p + acid_discrete_coord_flip()
        }
        ## Hide sample name legend.
        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(color = "none")
        }
        ## Return.
        p
    }



#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    f = "plot5Prime3PrimeBias",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plot5Prime3PrimeBias,bcbioRNASeq`
)
