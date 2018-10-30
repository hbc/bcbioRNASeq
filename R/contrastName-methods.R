#' Contrast Name
#'
#' @name contrastName
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `string`. Contrast name.
#'
#' @examples
#' data(deseq)
#'
#' ## DESeqAnalysis ====
#' contrastName(deseq, results = 1L)
#'
#' ## DESeqResults ====
#' object <- as(deseq, "DESeqResults")
#' contrastName(object)
NULL



contrastName.DESeqAnalysis <-  # nolint
    function(object, results) {
        do.call(
            what = contrastName,
            args = list(
                object = .matchResults(object, results)
            )
        )
    }



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



contrastName.DESeqResultsTables <-  # nolint
    function(object) {
        contrastName(slot(object, name = "results"))
    }



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqAnalysis"),
    definition = contrastName.DESeqAnalysis
)



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqResults"),
    definition = contrastName.DESeqResults
)



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqResultsTables"),
    definition = contrastName.DESeqResultsTables
)
