#' Design Formula Accessor
#'
#' @name design
#' @author Michael Steinbaugh
#'
#' @importFrom DESeq2 design design<-
#'
#' @inherit DESeq2::design
#'
#' @param object Object.
#'
#' @seealso
#' - [BiocGenerics::design].
#' - [DESeq2::design].
#'
#' @return `formula`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' design(bcb) %>% print(showEnv = FALSE)
#'
#' # Assignment support
#' design(bcb) <- ~group
#' design(bcb) %>% print(showEnv = FALSE)
NULL



# Methods ======================================================================
#' @rdname design
#' @export
setMethod(
    "design",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        dds <- assays(object)[["dds"]]
        assert_is_all_of(dds, "DESeqDataSet")
        validObject(dds)
        design(dds)
    }
)



#' @rdname design
#' @importFrom DESeq2 DESeq
#' @export
setMethod(
    "design<-",
    signature(
        object = "bcbioRNASeq",
        value = "formula"
    ),
    function(object, value) {
        dds <- assays(object)[["dds"]]
        assert_is_all_of(dds, "DESeqDataSet")
        design(dds) <- value
        dds <- DESeq(dds)
        validObject(dds)
        assays(object)[["dds"]] <- dds
        validObject(object)
        object
    }
)
