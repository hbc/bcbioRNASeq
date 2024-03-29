#' @name plotMeanSd
#' @author Michael Steinbaugh, Lorena Patano
#' @inherit AcidGenerics::plotMeanSd
#' @note Requires the vsn package to be installed.
#' @note Updated 2022-10-24.
#'
#' @inheritParams bcbioRNASeq
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @param lineColor `character(1)`.
#' Line color.
#'
#' @details
#' `vsn::meanSdPlot()` wrapper that plots count transformations on a log2 scale.
#'
#' - DESeq2 log2: log2 library size factor-adjusted normalized counts.
#' - DESeq2 rlog: **r**egularized **log** transformation.
#' - DESeq2 VST: **v**ariance-**s**tabilizing **t**ransformation.
#' - edgeR log2 TMM: log2 **t**rimmed **m**ean of **M**-values transformation.
#'
#' @seealso
#' - `vsn::meanSdPlot()`.
#' - `DESeq2::DESeq()`.
#' - `DESeq2::rlog()`.
#' - `DESeq2::varianceStabilizingTransformation()`.
#' - `edgeR::calcNormFactors()`.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' if (requireNamespace("vsn", quietly = TRUE)) {
#'     plotMeanSd(bcb)
#' }
NULL



## Updated 2022-03-07.
`plotMeanSd,bcbioRNASeq` <- # nolint
    function(object,
             fill = ggplot2::scale_fill_gradient(
                 low = AcidPlots::lightPalette[["gray"]],
                 high = AcidPlots::lightPalette[["purple"]]
             ),
             lineColor = AcidPlots::lightPalette[["orange"]],
             legend = getOption(x = "acid.legend", default = TRUE)) {
        validObject(object)
        assert(
            requireNamespaces("vsn"),
            isFlag(legend),
            isGgscale(fill, scale = "continuous", aes = "fill", nullOk = TRUE),
            isString(lineColor, nullOk = TRUE)
        )
        ## Determine which genes are non-zero, and should be included in plot.
        raw <- counts(object, normalized = FALSE)
        nonzero <- rowSums(raw) > 0L
        args <- list(
            "sf" = list(log2 = FALSE, title = "sf"),
            "rlog" = list(log2 = TRUE, title = "rlog"),
            "vst" = list(log2 = TRUE, title = "vst"),
            "tmm" = list(log2 = FALSE, title = "tmm"),
            "rle" = list(log2 = FALSE, title = "rle")
        )
        normalized <- names(args)
        log2 <- vapply(
            X = args,
            FUN = function(x) x[["log2"]],
            FUN.VALUE = logical(1L)
        )
        ## Get the requested counts from object.
        assays <- Map(
            normalized = normalized,
            log2 = log2,
            MoreArgs = list(
                "object" = object,
                "nonzero" = nonzero
            ),
            f = function(normalized, log2, object, nonzero) {
                suppressMessages({
                    mat <- tryCatch(
                        expr = counts(object, normalized = normalized),
                        error = function(e) NULL
                    )
                })
                if (!is.matrix(mat)) {
                    return(NULL)
                }
                assert(identical(length(nonzero), nrow(mat)))
                mat <- mat[nonzero, , drop = FALSE]
                if (!isTRUE(log2)) {
                    mat <- log2(mat + 1L)
                }
                mat
            }
        )
        assays <- Filter(f = Negate(is.null), x = assays)
        titles <- vapply(
            X = args,
            FUN = function(x) x[["title"]],
            FUN.VALUE = character(1L)
        )
        titles <- titles[names(assays)]
        plotlist <- Map(
            assay = assays,
            title = titles,
            MoreArgs = list(
                "fill" = fill,
                "legend" = legend,
                "lineColor" = lineColor
            ),
            f = function(assay, title, fill, legend, lineColor) {
                p <- vsn::meanSdPlot(
                    x = assay,
                    ranks = TRUE,
                    plot = FALSE,
                    bins = 50L
                )
                p <- p[["gg"]] + ggtitle(paste(title, "(log2)"))
                ## Improve the fill aesthetics.
                if (is(fill, "ScaleContinuous")) {
                    suppressMessages(p <- p + fill)
                }
                ## Improve the line aesthetics.
                p[["layers"]][[2L]][["aes_params"]][["colour"]] <- lineColor
                p[["layers"]][[2L]][["aes_params"]][["size"]] <- 1L
                p
            }
        )
        ## Remove the plot (color) legend, if desired.
        if (!isTRUE(legend)) {
            plotlist <- lapply(plotlist, function(p) {
                p <- p + theme(legend.position = "none")
            })
        }
        wrap_plots(plotlist)
    }



#' @rdname plotMeanSd
#' @export
setMethod(
    f = "plotMeanSd",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotMeanSd,bcbioRNASeq`
)
