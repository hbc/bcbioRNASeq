#' Contrast Name
#'
#' @name contrastName
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return Contrast name `string`.
#'
#' @examples
#' load(system.file("extdata/res.rda", package = "bcbioRNASeq"))
#'
#' # DESeqResults
#' contrastName(res)
NULL



# Methods ======================================================================
#' @rdname contrastName
#' @export
setMethod(
    "contrastName",
    signature("DESeqResults"),
    function(object) {
        validObject(object)
        contrast <- mcols(object)[2L, 2L]
        assert_is_character(contrast)
        contrast %>%
            gsub("^.*:\\s", "", .) %>%
            gsub("_", " ", .) %>%
            # Improve appearance for difference of differences
            gsub("\\+", " \\+\n    ", .)
    }
)
