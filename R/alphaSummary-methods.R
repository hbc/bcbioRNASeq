#' Alpha Level Cutoff Summary Statistics
#'
#' Quickly generate a summary table of various alpha level cutoffs, for use in
#' an R Markdown report.
#'
#' @note Use either `contrast` or `name` to specify the desired contrast.
#'
#' @name alphaSummary
#' @family Differential Expression Functions
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams general
#' @inheritParams DESeq2::results
#' @param alpha `numeric`. Multiple alpha cutoffs.
#' @param caption `string` or `NULL`. Table caption. If set `NULL`, will be
#'   generated automatically from the contrast used.
#'
#' @return `kable`.
#'
#' @seealso
#' - [DESeq2::results()].
#' - [DESeq2::resultsNames()].
#'
#' @examples
#' # DESeqDataSet ====
#' object <- deseq_small@data
#' design(object)
#' resultsNames(object)
#' alphaSummary(object, contrast = c("treatment", "folic_acid", "control"))
#' alphaSummary(object, name = "treatment_folic_acid_vs_control")
NULL



#' @rdname alphaSummary
#' @export
setMethod(
    f = "alphaSummary",
    signature = signature("DESeqDataSet"),
    definition = function(
        object,
        alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
        contrast,
        name,
        caption = NULL
    ) {
        validObject(object)
        assert_is_numeric(alpha)
        assertIsAStringOrNULL(caption)
        # Require either `contrast` or `name`.
        if (
            (missing(contrast) && missing(name)) ||
            (!missing(contrast) && !missing(name))
        ) {
            stop(
                "Specify either `contrast` or `name` (but not both)",
                call. = FALSE
            )
        } else if (!missing(contrast)) {
            assert_is_character(contrast)
            name <- NULL
        } else if (!missing(name)) {
            assert_is_a_string(name)
            contrast <- NULL
        }

        # Generate an automatic caption.
        if (is.null(caption)) {
            if (!is.null(contrast)) {
                caption <- paste(contrast, collapse = " ")
            } else if (!is.null(name)) {
                caption <- name
            }
        }

        data <- sapply(
            X = alpha,
            FUN = function(alpha) {
                args <- list(
                    object = object,
                    contrast = contrast,
                    name = name,
                    alpha = alpha
                )
                args <- Filter(Negate(is.null), args)
                results <- do.call(what = results, args = args)
                output <- capture.output(summary(results))
                # Subset the lines of interest from summary.
                # Keep only the summary lines that contain a colon.
                output <- output[grepl(" : ", output)]
                # Extract the values after the colon in summary.

                x <- str_match(
                    string = output,
                    pattern = "^(.+)\\s\\:\\s(.+)$"
                )
                names <- gsub("\\s+$", "", x[, 2L])
                values <- x[, 3L]
                names(values) <- names
                values
            }
        )
        # Coerce to data.frame.
        data <- data %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            set_colnames(as.character(alpha))

        kable(data, caption = caption)

    }
)
