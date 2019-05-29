#' @name contrastName
#' @inherit bioverbs::contrastName
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `string`. Contrast name.
#'
#' @examples
#' # DESeqResults ====
#' contrastName(res_small)
NULL



#' @rdname contrastName
#' @name contrastName
#' @importFrom bioverbs contrastName
#' @usage contrastName(object, ...)
#' @export
NULL



contrastName.DESeqResults <-  # nolint
    function(object) {
        validObject(object)
        contrast <- mcols(object)[2L, "description"]
        assert_is_character(contrast)
        contrast %>%
            gsub("^.*:\\s", "", .) %>%
            gsub("_", " ", .) %>%
            # Improve appearance for difference of differences
            gsub("\\+", " \\+\n    ", .)
    }



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqResults"),
    definition = contrastName.DESeqResults
)
