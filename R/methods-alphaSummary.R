#' Alpha Level Cutoff Summary Statistics
#'
#' @note `bcbioRNASeq` class does currently support contrast definitions, since
#'   the object contains an internal `DESeqDataSet` with an empty design
#'   formula.
#'
#' @name alphaSummary
#' @family Differential Expression Functions
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams general
#' @inheritParams DESeq2::results
#' @param alpha Numeric vector of multiple alpha cutoffs.
#' @param caption *Optional.* Character string to use as a caption.
#' @param ... *Optional.* Passthrough arguments to [DESeq2::results()]. Use
#'   either `contrast` or `name` arguments to define the desired contrast.
#'
#' @return `kable`.
#'
#' @seealso [DESeq2::results()].
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds_small.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq ====
#' alphaSummary(bcb_small)
#'
#' # DESeqDataSet ====
#' resultsNames(dds_small)
#' alphaSummary(dds_small, contrast = c("treatment", "folic_acid", "control"))
#' alphaSummary(dds_small, name = "treatment_folic_acid_vs_control")
NULL



# Methods ======================================================================
#' @rdname alphaSummary
#' @export
setMethod(
    "alphaSummary",
    signature("DESeqDataSet"),
    function(
        object,
        alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
        caption = NULL,
        ...
    ) {
        validObject(object)
        assert_is_numeric(alpha)
        assertIsAStringOrNULL(caption)
        dots <- list(...)

        # Abort on empty design formula
        if (design(object) == ~ 1) {  # nolint
            warn("Empty design formula detected")
        }

        # Generate an automatic caption
        if (is.null(caption)) {
            if (!is.null(dots[["contrast"]])) {
                caption <- paste(dots[["contrast"]], collapse = " ")
            } else if (!is.null(dots[["name"]])) {
                caption <- dots[["name"]]
            }
        }

        dflist <- lapply(alpha, function(x) {
            output <- capture.output(
                summary(results(object, ..., alpha = x))
            )
            # Subset the lines of interest from summary
            output <- output[4L:8L]
            # Extract the values after the colon in summary
            values <- vapply(
                X = output,
                FUN = function(x) {
                    gsub("^.+\\:\\s(.+)\\s$", "\\1", x)
                },
                FUN.VALUE = "character"
            )
            data.frame(alpha = values)
        })

        dflist %>%
            bind_cols() %>%
            set_colnames(alpha) %>%
            set_rownames(c(
                "LFC > 0 (up)",
                "LFC < 0 (down)",
                "outliers",
                "low counts",
                "cutoff"
            )) %>%
            kable(caption = caption)
    }
)



#' @rdname alphaSummary
#' @export
setMethod(
    "alphaSummary",
    signature("bcbioRNASeq"),
    function(
        object,
        alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
        caption = NULL,
        ...
    ) {
        validObject(object)
        dds <- as(object, "DESeqDataSet")
        alphaSummary(
            object = dds,
            alpha = alpha,
            caption = caption,
            ...
        )
    }
)
