#' Plot Sexually Dimorphic Gender Markers
#'
#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams plotGene
#' @inheritParams plotTotalReads
#'
#' @param organism *Optional*. Organism name. Should be detected automatically,
#'   unless a spike-in FASTA sequence is provided containing a gene identifier
#'   that is first alphabetically in the count matrix rownames.
#'
#' @return [ggplot].
#'
#' @examples
#' # bcbioRNASeq
#' # Use F1000 workflow example dataset
#' # The minimal example inside the package doesn't have dimorphic genes
#' loadRemoteData(
#'     file.path(
#'         "http://bcbiornaseq.seq.cloud",
#'         "f1000v1",
#'         "data",
#'         "bcb.rda"),
#'     quiet = TRUE)
#'
#' # bcbioRNASeq
#' plotGenderMarkers(bcb)
#' plotGenderMarkers(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     color = NULL)
#'
#' # DESeqDataSet
#' dds <- bcbio(bcb, "DESeqDataSet")
#' plotGenderMarkers(dds)
NULL



# Constructors =================================================================
#' @importFrom dplyr filter left_join pull
#' @importFrom ggplot2 aes_string expand_limits geom_jitter ggplot labs
#' @importFrom stats setNames
#' @importFrom tibble rownames_to_column
#' @importFrom viridis scale_color_viridis
.plotGenderMarkers <- function(
    object,
    interestingGroups = "sampleName",
    organism,
    metadata,
    countsAxisLabel = "counts",
    color = viridis::scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    if (isTRUE(title)) {
        title <- "gender markers"
    } else if (!is.character(title)) {
        title <- NULL
    }

    # Load the relevant internal gender markers data
    markers <- get("genderMarkers", envir = loadNamespace("bcbioRNASeq"))
    if (!camel(organism) %in% names(markers)) {
        warn(paste(
            "Organism",
            paste0("(", organism, ")"),
            "is not supported"
        ))
        return(invisible(NULL))
    }
    # Convert the organism name from full latin to shorthand (e.g. hsapiens)
    markers <- markers[[camel(organism)]]

    # Ensembl identifiers
    ensgene <- markers %>%
        filter(.data[["include"]] == TRUE) %>%
        pull("ensgene") %>%
        sort() %>%
        unique()

    if (!all(ensgene %in% rownames(object))) {
        warn("Missing gender markers in count matrix")
        return(invisible(NULL))
    }

    data <- object %>%
        .[ensgene, , drop = FALSE] %>%
        # This will coerce rownames to a column named `rowname`. We will rename
        # this to `ensgene` after melting the counts.
        as.data.frame() %>%
        rownames_to_column() %>%
        # For `melt()`, can also declare `measure.vars` here instead of using
        # `setNames()`. If you don't set `id`, function will output a message.
        melt(id = 1L) %>%
        setNames(c("ensgene", "sampleName", "counts")) %>%
        left_join(markers, by = "ensgene") %>%
        left_join(metadata, by  = "sampleName") %>%
        uniteInterestingGroups(interestingGroups)

    p <- ggplot(
        data,
        mapping = aes_string(
            x = "symbol",
            y = "counts",
            color = "interestingGroups",
            shape = "chromosome")
    ) +
        geom_jitter(size = 4L) +
        expand_limits(y = 0L) +
        labs(
            title = title,
            x = "gene",
            y = countsAxisLabel,
            color = paste(interestingGroups, collapse = ":\n"))

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    p
}



# Methods ======================================================================
#' @rdname plotGenderMarkers
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenderMarkers",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        color = viridis::scale_color_viridis(discrete = TRUE),
        title = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotGenderMarkers(
            object = tpm(object),
            interestingGroups = interestingGroups,
            organism = metadata(object)[["organism"]],
            metadata = sampleMetadata(object),
            countsAxisLabel = "transcripts per million (tpm)",
            color = color,
            title = title)
    })



#' @rdname plotGenderMarkers
#' @importFrom bcbioBase detectOrganism
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqDataSet"),
    function(
        object,
        interestingGroups = "sampleName",
        organism,
        color = viridis::scale_color_viridis(discrete = TRUE),
        title = TRUE) {
        counts <- counts(object, normalized = TRUE)
        if (missing(organism)) {
            organism <- rownames(counts) %>%
                .[[1L]] %>%
                detectOrganism()
        }
        .plotGenderMarkers(
            counts,
            interestingGroups = interestingGroups,
            organism = organism,
            metadata = sampleMetadata(object),
            countsAxisLabel = "normalized counts",
            color = color,
            title = title)
    })



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("matrix"),
    .plotGenderMarkers)
