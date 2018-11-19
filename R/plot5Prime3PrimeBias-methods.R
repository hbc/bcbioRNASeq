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



# bcbioRNASeq ==================================================================
plot5Prime3PrimeBias.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0L,
        fill,
        flip,
        title = "5'->3' bias"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assert_is_a_number(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsStringOrNULL(title)

        data <- metrics(object)

        # The formatting of this column can vary depending on the version of
        # `camel()` used. This change was added in v0.2.7.
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
                y = "5'->3' bias",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            p <- p + basejump_geom_abline(yintercept = limit)  # nocov
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

formals(plot5Prime3PrimeBias.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plot5Prime3PrimeBias.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    f = "plot5Prime3PrimeBias",
    signature = signature("bcbioRNASeq"),
    definition = plot5Prime3PrimeBias.bcbioRNASeq
)
