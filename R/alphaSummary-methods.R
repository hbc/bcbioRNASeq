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
#' @export
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
        contrast = NULL,
        name = NULL,
        caption = NULL
    ) {
        validObject(object)
        assert_is_numeric(alpha)
        assert_is_any_of(contrast, c("character", "NULL"))
        assertIsAStringOrNULL(name)
        assertIsAStringOrNULL(caption)

        # Either `contrast` or `name`.
        # If neither are defined, we're checking the intercept.
        if (!is.null(contrast) && !is.null(name)) {
            stop(
                "Specify either `contrast` or `name` (but not both)",
                call. = FALSE
            )
        }

        # Generate an automatic caption.
        if (is.null(caption)) {
            if (!is.null(contrast)) {
                caption <- paste(contrast, collapse = " ")
            } else if (!is.null(name)) {
                caption <- name
            } else {
                caption <- resultsNames(object)[[1L]]
            }
        }
        message(paste(
            caption,
            paste("Alpha cutoffs:", toString(alpha)),
            sep = "\n"
        ))

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
