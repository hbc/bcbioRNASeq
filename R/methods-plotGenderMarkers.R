#' Plot Sexually Dimorphic Gender Markers
#'
#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotGene
#' @inheritParams plotTotalReads
#'
#' @param organism *Optional.* Organism name. Should be detected automatically,
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
#' @importFrom basejump camel
#' @importFrom dplyr filter left_join pull
#' @importFrom ggplot2 aes_string expand_limits geom_jitter ggplot labs
#' @importFrom magrittr set_colnames
#' @importFrom tibble rownames_to_column
.plotGenderMarkers <- function(
    object,
    interestingGroups = "sampleName",
    organism,
    metadata,
    countsAxisLabel = "counts",
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_matrix(object)
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_a_string(organism)
    assert_is_data.frame(metadata)
    assert_is_a_string(countsAxisLabel)
    .assert_formal_discrete_scale(color)
    .assert_formal_title(title)
    
    if (isTRUE(title)) {
        title <- "gender markers"
    } else if (!is.character(title)) {
        title <- NULL
    }

    # Load the relevant internal gender markers data
    markers <- get("genderMarkers", envir = loadNamespace("bcbioRNASeq"))
    assert_is_subset(camel(organism), names(markers))
    markers <- markers[[camel(organism)]]

    # Ensembl identifiers
    ensgene <- markers %>%
        .[.[["include"]] == TRUE, "ensgene", drop = TRUE] %>%
        sort() %>%
        unique()
    assert_is_subset(ensgene, rownames(object))

    data <- object %>%
        .[ensgene, , drop = FALSE] %>%
        # This will coerce rownames to a column named `rowname`. We will rename
        # this to `ensgene` after melting the counts.
        as.data.frame() %>%
        rownames_to_column() %>%
        # For `melt()` can also declare `measure.vars` here instead
        melt(id = 1L) %>%
        set_colnames(c("ensgene", "sampleName", "counts")) %>%
        left_join(markers, by = "ensgene") %>%
        left_join(metadata, by  = "sampleName") %>%
        uniteInterestingGroups(interestingGroups)

    p <- ggplot(
        data = data,
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



.plotGenderMarkers.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    normalized = "tpm",
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_a_string(normalized)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotGenderMarkers(
        object = counts(object, normalized = normalized),
        interestingGroups = interestingGroups,
        organism = metadata(object)[["organism"]],
        metadata = sampleMetadata(object),
        countsAxisLabel = normalized,
        color = color,
        title = title)
}



.plotGenderMarkers.DESeqDataSet <- function(  # nolint
    object,
    interestingGroups = "sampleName",
    organism,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    # Passthrough: interestingGroups, color, title
    counts <- counts(object, normalized = TRUE)
    if (missing(organism)) {
        organism <- detectOrganism(counts)
    }
    .plotGenderMarkers(
        object = counts,
        interestingGroups = interestingGroups,
        organism = organism,
        metadata = sampleMetadata(object),
        countsAxisLabel = "normalized counts",
        color = color,
        title = title)
}



# Methods ======================================================================
#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("bcbioRNASeq"),
    .plotGenderMarkers.bcbioRNASeq)



#' @rdname plotGenderMarkers
#' @importFrom basejump detectOrganism
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqDataSet"),
    .plotGenderMarkers.DESeqDataSet)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("matrix"),
    .plotGenderMarkers)
