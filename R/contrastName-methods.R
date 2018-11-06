#' @name contrastName
#' @inherit basejump::contrastName
#' @author Michael Steinbaugh
#'
#' @inheritParams params
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



#' @importFrom basejump contrastName
#' @aliases NULL
#' @export
basejump::contrastName



# DESeqAnalysis ================================================================
contrastName.DESeqAnalysis <-  # nolint
    function(object, results) {
        do.call(
            what = contrastName,
            args = list(
                object = .matchResults(object, results)
            )
        )
    }



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqAnalysis"),
    definition = contrastName.DESeqAnalysis
)



# DESeqResults =================================================================
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



# DESeqResultsTables ===========================================================
contrastName.DESeqResultsTables <-  # nolint
    function(object) {
        contrastName(slot(object, name = "results"))
    }



#' @rdname contrastName
#' @export
setMethod(
    f = "contrastName",
    signature = signature("DESeqResultsTables"),
    definition = contrastName.DESeqResultsTables
)
