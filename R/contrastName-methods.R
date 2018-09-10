#' Contrast Name
#'
#' @name contrastName
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `string`. Contrast name.
#'
#' @examples
#' # DESeqAnalysis ====
#' contrastName(deseq_small, results = 1L)
#'
#' # DESeqResults ====
#' object <- deseq_small@results[[1L]]
#' contrastName(object)
NULL



#' @rdname contrastName
#' @export
setMethod(
    "contrastName",
    signature("DESeqAnalysis"),
    function(object, results) {
        do.call(
            what = contrastName,
            args = list(
                object = .matchResults(object, results)
            )
        )
    }
)



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



#' @rdname contrastName
#' @export
setMethod(
    "contrastName",
    signature("DESeqResultsTables"),
    function(object) {
        contrastName(object@all)
    }
)
