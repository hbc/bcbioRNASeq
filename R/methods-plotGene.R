#' Plot Individual Genes
#'
#' @rdname plotGene
#' @name plotGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotGene
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param format Ensembl identifier format. Supports `symbol` (**default**) or
#'   `ensgene`.
#'
#' @return [ggplot].
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' data(bcb)
#' genes <- c("Sulf1", "Phf3")
#'
#' # bcbioRNASeq
#' plotGene(bcb, gene = genes)
#'
#' \dontrun{
#' plotGene(bcb, gene = genes, interestingGroups = "group")
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_color_viridis
.plotGene <- function(
    counts,
    gene,
    metadata,
    interestingGroups = "sampleName") {
    metadata <- as.data.frame(metadata)
    sapply(seq_along(gene), function(a) {
        ensgene <- gene[[a]]
        symbol <- names(gene)[[a]]
        df <- data.frame(
            x = colnames(counts),
            y = counts[ensgene, ],
            color = metadata[[interestingGroups]])
        p <- ggplot(
            df,
            mapping = aes_string(
                x = "x",
                y = "y",
                color = "color")
        ) +
            geom_point(size = 4) +
            theme(
                axis.text.x = element_text(angle = 90)) +
            labs(title = symbol,
                 x = "sample",
                 y = "counts",
                 color = interestingGroups) +
            expand_limits(y = 0) +
            scale_color_viridis(discrete = TRUE)
        show(p)
    }) %>%
        invisible()
}



# Methods ====
#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        normalized = "tpm",
        gene,
        format = "symbol") {
        if (!format %in% c("ensgene", "symbol")) {
            stop("Unsupported gene identifier format", call. = FALSE)
        }
        if (is.logical(normalized)) {
            warning(paste(
                "Explicit format of normalized counts format is recommended"
            ), call. = FALSE)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }

        counts <- counts(object, normalized = normalized)
        metadata <- colData(object)

        # Match unique gene identifier with name (gene symbol) using the
        # internally stored Ensembl annotations saved in the run object
        gene2symbol <- metadata(object) %>%
            .[["annotable"]] %>%
            .[.[[format]] %in% gene, , drop = FALSE] %>%
            .[, c("symbol", "ensgene")]
        gene <- gene2symbol[["ensgene"]]
        names(gene) <- gene2symbol[["symbol"]]

        .plotGene(
            counts = counts,
            gene = gene,
            metadata = metadata,
            interestingGroups = interestingGroups)
    })
