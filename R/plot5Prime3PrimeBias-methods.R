#' @name plot5Prime3PrimeBias
#' @author Michael Steinbaugh
#' @inherit bioverbs::plot5Prime3PrimeBias
#' @note Updated 2019-09-13.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plot5Prime3PrimeBias(bcb)
NULL



#' @rdname plot5Prime3PrimeBias
#' @name plot5Prime3PrimeBias
#' @importFrom bioverbs plot5Prime3PrimeBias
#' @usage plot5Prime3PrimeBias(object, ...)
#' @export
NULL



## FIXME Rework labels

## Updated 2019-09-13.
`plot5Prime3PrimeBias,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        color,
        labels = list(
            title = "5'->3' bias",
            subtitle = NULL,
            x = NULL,
            y = "5'->3' bias"
        ),
        flip
    ) {
        validObject(object)
        assert(
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE),
            is.list(labels),
            areSetEqual(
                x = names(labels),
                y = names(eval(formals()[["labels"]]))
            ),
            isFlag(flip)
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        ## The formatting of this column can vary depending on the version of
        ## `camelCase()` used. This grep match fix was added in v0.2.7.
        yCol <- grep(
            pattern = ".+5.+3bias$",
            x = colnames(data),
            ignore.case = TRUE,
            value = TRUE
        )
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym(yCol),
                color = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 3L) +
            acid_geom_abline(yintercept = 1L)
        ## Labels.
        if (is.list(labels)) {
            labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
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
