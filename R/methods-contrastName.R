#' Contrast Name
#'
#' @name contrastName
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return Contrast name `string`.
#'
#' @examples
#' load(system.file("extdata/res_small.rda", package = "bcbioRNASeq"))
#'
#' # DESeqResults ====
#' contrastName(res_small)
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
