#' @name plot5Prime3PrimeBias
#' @author Michael Steinbaugh
#' @inherit bioverbs::plot5Prime3PrimeBias
#' @note Updated 2019-08-07.
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



## Updated 2019-07-23.
`plot5Prime3PrimeBias,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        color,
        flip,
        title = "5'->3' bias"
    ) {
        validObject(object)
        assert(
            isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE),
            isFlag(flip),
            isString(title, nullOK = TRUE)
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
                colour = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 3L) +
            labs(
                title = title,
                x = NULL,
                y = "5'->3' bias",
                colour = paste(interestingGroups, collapse = ":\n")
            ) +
            acid_geom_abline(yintercept = 1L)

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        if (isTRUE(flip)) {
            p <- acid_coord_flip(p)
        }

        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(color = FALSE)
        }

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
