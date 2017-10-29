#' Plot Sexually Dimorphic Gender Markers
#'
#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams plotGene
#'
#' @param organism *Optional*. Organism name. Should be detected automatically,
#'   unless a spike-in FASTA sequence is provided containing a gene identifier
#'   that is first alphabetically in the count matrix rownames.
#' @param ylab Y-axis label.
#'
#' @return [ggplot].
#'
#' @examples
#' # Use F1000 workflow example dataset
#' # The minimal example inside the package doesn't have dimorphic genes
#' \dontrun{
#' file.path(
#'     "https://github.com",
#'     "hbc",
#'     "bcbioRNASeq",
#'     "raw",
#'     "f1000v1",
#'     "data",
#'     "bcb.rda") %>%
#'     loadRemoteData()
#' plotGenderMarkers(bcb)
#' }
NULL



# Constructors ====
#' @importFrom dplyr left_join pull
#' @importFrom ggplot2 aes_string expand_limits geom_jitter ggplot labs
#' @importFrom stats setNames
#' @importFrom viridis scale_color_viridis
.plotGenderMarkers <- function(
    object,
    organism,
    ylab = "counts",
    color = scale_color_viridis(discrete = TRUE)) {
    # Load the relevant internal gender markers data
    envir <- loadNamespace("bcbioRNASeq")
    if (organism == "Homo sapiens") {
        markers <- get("genderMarkersHsapiens", envir = envir)
    } else if (organism == "Mus musculus") {
        markers <- get("genderMarkersMmusculus", envir = envir)
    } else {
        warning("Unsupported organism", call. = FALSE)
        return(NULL)
    }

    # Ensembl identifiers
    ensgene <- markers %>%
        .[.[["include"]] == TRUE, , drop = FALSE] %>%
        pull("ensgene") %>%
        sort() %>%
        unique()

    if (!all(ensgene %in% rownames(object))) {
        warning("Missing gender markers in count matrix", call. = FALSE)
        return(NULL)
    }

    counts <- object %>%
        .[ensgene, , drop = FALSE] %>%
        # This will coerce rownames to a column named `rowname`. We will rename
        # this to `ensgene` after melting the counts.
        as("tibble") %>%
        # For `melt()`, can also declare `measure.vars` here instead of using
        # `setNames()`. If you don't set `id`, function will output a message.
        melt(id = 1) %>%
        setNames(c("ensgene", "sampleName", "counts")) %>%
        left_join(markers, by = "ensgene")
    p <- ggplot(
        counts,
        mapping = aes_string(
            x = "symbol",
            y = "counts",
            color = "sampleName",
            shape = "chromosome")
    ) +
        geom_jitter(size = 4) +
        expand_limits(y = 0) +
        labs(title = "gender markers",
             x = "gene",
             y = ylab)
    if (!is.null(color)) {
        p <- p + color
    }
    p
}



# Methods ====
#' @rdname plotGenderMarkers
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenderMarkers",
    signature("bcbioRNASeq"),
    function(
        object,
        color = scale_color_viridis(discrete = TRUE)) {
        counts <- tpm(object)
        organism <- metadata(object)[["organism"]]
        ylab <- "transcripts per million (tpm)"
        .plotGenderMarkers(
            counts,
            organism = organism,
            ylab = ylab,
            color = color)
    })



#' @rdname plotGenderMarkers
#' @importFrom basejump detectOrganism
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqDataSet"),
    function(
        object,
        organism,
        color = scale_color_viridis(discrete = TRUE)) {
        counts <- counts(object, normalized = TRUE)
        if (missing(organism)) {
            organism <- rownames(counts) %>%
                .[[1]] %>%
                detectOrganism()
        }
        ylab <- "normalized counts"
        .plotGenderMarkers(
            counts,
            organism = organism,
            ylab = ylab,
            color = color)
    })



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("matrix"),
    .plotGenderMarkers)
