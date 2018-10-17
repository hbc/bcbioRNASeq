#' Contrast Name
#'
#' @name contrastName
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `string`. Contrast name.
#'
#' @examples
#' data(deseq_small)
#' 
#' ## DESeqAnalysis ====
#' contrastName(deseq_small, results = 1L)
#'
#' ## DESeqResults ====
#' object <- deseq_small@results[[1L]]
#' contrastName(object)
NULL



.contrastName.DESeqAnalysis <-  # nolint
    function(object, results) {
        do.call(
            what = contrastName,
            args = list(
                object = .matchResults(object, results)
            )
        )
    }



.contrastName.DESeqResults <-  # nolint
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



.contrastName.DESeqResultsTables <-  # nolint
    function(object) {
        contrastName(slot(object, name = "results"))
    }



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqAnalysis"),
    definition = .contrastName.DESeqAnalysis
)



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqResults"),
    definition = .contrastName.DESeqResults
)



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqResultsTables"),
    definition = .contrastName.DESeqResultsTables
)
