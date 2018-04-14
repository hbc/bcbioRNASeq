#' Show an Object
#'
#' @name show
#' @family S4 Class Definition
#' @author Michael Steinbuagh
#'
#' @inherit methods::show
#'
#' @examples
#' show(bcb_small)
NULL



# Methods ======================================================================
#' @rdname show
#' @export
setMethod(
    "show",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)

        cat(
            paste(class(object), metadata(object)[["version"]]),
            paste("Samples:", ncol(object)),
            paste0(
                sub(
                    pattern = "^([a-z])",
                    replacement = "\\U\\1",
                    x = metadata(object)[["level"]],
                    perl = TRUE
                ),
                ": ",
                nrow(object)
            ),
            paste("Organism:", metadata(object)[["organism"]]),
            sep = "\n"
        )

        # rowRanges
        m <- metadata(object)[["rowRangesMetadata"]]
        if (is.data.frame(m) && length(m)) {
            cat(
                paste(
                    "AnnotationHub:",
                    m[m[["name"]] == "id", "value", drop = TRUE]
                ),
                paste(
                    "Ensembl Release:",
                    m[m[["name"]] == "ensembl_version", "value", drop = TRUE]
                ),
                paste(
                    "Genome Build:",
                    m[m[["name"]] == "genome_build", "value", drop = TRUE]
                ),
                sep = "\n"
            )
        }

        cat(
            paste("Upload Dir:", metadata(object)[["uploadDir"]]),
            paste("Upload Date:", metadata(object)[["runDate"]]),
            paste("R Load Date:", metadata(object)[["date"]]),
            sep = "\n"
        )
    }
)
