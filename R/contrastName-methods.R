#' Contrast Name
#'
#' @name contrastName
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `string`. Contrast name.
#'
#' @examples
#' # DESeqResults ====
#' contrastName(res_small)
NULL



#' @rdname contrastName
#' @export
setMethod(
    "contrastName",
    signature("DESeqResults"),
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
)
