#' @name plot5Prime3PrimeBias
#' @inherit basejump::plot5Prime3PrimeBias
#' @author Michael Steinbaugh
#'
#' @inheritParams params
#' @inheritParams basejump::params
#'
#' @examples
#' data(bcb)
#' plot5Prime3PrimeBias(bcb)
NULL



#' @importFrom basejump plot5Prime3PrimeBias
#' @aliases NULL
#' @export
basejump::plot5Prime3PrimeBias



plot5Prime3PrimeBias.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        color,
        flip,
        title = "5'->3' bias"
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        # TODO Update isGGScale and isString assertions to optionally allow
        # NULL for convenience. This should be FALSE by default.
        assert(
            isGGScale(color, scale = "discrete", aes = "colour") ||
                is.null(color),
            isFlag(flip),
            isString(title) || is.null(title)
        )

        data <- metrics(object)

        # The formatting of this column can vary depending on the version of
        # `camel` used. This grep match fix was added in v0.2.7.
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
            basejump_geom_abline(yintercept = 1L)

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        if (isTRUE(flip)) {
            p <- p + coord_flip()
        }

        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(color = FALSE)
        }

        p
    }

formals(plot5Prime3PrimeBias.bcbioRNASeq)[["color"]] <-
    formalsList[["color.discrete"]]
formals(plot5Prime3PrimeBias.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    f = "plot5Prime3PrimeBias",
    signature = signature("bcbioRNASeq"),
    definition = plot5Prime3PrimeBias.bcbioRNASeq
)
