#' @name plot5Prime3PrimeBias
#' @author Michael Steinbaugh
#' @inherit AcidGenerics::plot5Prime3PrimeBias
#' @note Updated 2019-09-16.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plot5Prime3PrimeBias(bcb)
NULL



#' @rdname plot5Prime3PrimeBias
#' @name plot5Prime3PrimeBias
#' @importFrom AcidGenerics plot5Prime3PrimeBias
#' @usage plot5Prime3PrimeBias(object, ...)
#' @export
NULL



## Updated 2019-09-16.
`plot5Prime3PrimeBias,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        color,
        labels = list(
            title = "5'->3' bias",
            subtitle = NULL,
            sampleAxis = NULL,
            metricAxis = "5'->3' bias"
        ),
        flip
    ) {
        validObject(object)
        assert(
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE),
            isFlag(flip)
        )
        labels <- matchLabels(
            labels = labels,
            choices = eval(formals()[["labels"]])
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("x5x3Bias"),
                color = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 3L) +
            acid_geom_abline(yintercept = 1L)
        ## Labels.
        if (is.list(labels)) {
            labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
            names(labels)[names(labels) == "sampleAxis"] <- "x"
            names(labels)[names(labels) == "metricAxis"] <- "y"
            p <- p + do.call(what = labs, args = labels)
        }
        ## Color.
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
        ## Flip.
        if (isTRUE(flip)) {
            p <- acid_coord_flip(p)
        }
        ## Hide sample name legend.
        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(color = FALSE)
        }
        ## Return.
        p
    }

formals(`plot5Prime3PrimeBias,bcbioRNASeq`)[c("color", "flip")] <-
    formalsList[c("color.discrete", "flip")]



#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    f = "plot5Prime3PrimeBias",
    signature = signature("bcbioRNASeq"),
    definition = `plot5Prime3PrimeBias,bcbioRNASeq`
)
