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
#' dir <- "http://bcbiornaseq.seq.cloud/f1000v1/data"
#' loadRemoteData(
#'     c(
#'         file.path(dir, "bcb.rda"),
#'         file.path(dir, "rld.rda"),
#'         file.path(dir, "vst.rda")
#'     ),
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
#' plotGenderMarkers(dds, interestingGroups = "group")
#'
#' # DESeqTransform
#' plotGenderMarkers(rld, interestingGroups = "group")
#' plotGenderMarkers(vst, interestingGroups = "group")
NULL



# Constructors =================================================================
#' @importFrom basejump camel
#' @importFrom dplyr filter left_join pull
#' @importFrom ggplot2 aes_string expand_limits geom_jitter ggplot labs
#' @importFrom magrittr set_colnames
#' @importFrom tibble rownames_to_column
.plotGenderMarkers <- function(
    object,
    metadata,
    interestingGroups = "sampleName",
    organism,
    countsAxisLabel = "counts",
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_matrix(object)
    assert_is_data.frame(metadata)
    assert_formal_interesting_groups(metadata, interestingGroups)
    assert_is_a_string(organism)
    assert_is_a_string(countsAxisLabel)
    .assert_formal_scale_discrete(color)
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
    normalized = "rlog",
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



#' @importFrom basejump detectOrganism
.plotGenderMarkers.DESeqDataSet <- function(  # nolint
    object,
    interestingGroups = "sampleName",
    organism,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    # Passthrough: interestingGroups, color, title
    counts <- log2(counts(object, normalized = TRUE) + 1L)
    if (missing(organism)) {
        organism <- detectOrganism(counts)
    }
    .plotGenderMarkers(
        object = counts,
        interestingGroups = interestingGroups,
        organism = organism,
        metadata = sampleMetadata(object),
        countsAxisLabel = "log2 normalized counts",
        color = color,
        title = title)
}



#' @importFrom basejump detectOrganism
.plotGenderMarkers.DESeqTransform <- function(  # nolint
    object,
    interestingGroups = "sampleName",
    organism,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    # Passthrough: interestingGroups, color, title
    counts <- assay(object)
    if (missing(organism)) {
        organism <- detectOrganism(counts)
    }
    if ("rlogIntercept" %in% colnames(mcols(object))) {
        countsAxisLabel <- "rlog"
    } else {
        countsAxisLabel <- "vst"
    }
    .plotGenderMarkers(
        object = counts,
        interestingGroups = interestingGroups,
        organism = organism,
        metadata = sampleMetadata(object),
        countsAxisLabel = countsAxisLabel,
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
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqDataSet"),
    .plotGenderMarkers.DESeqDataSet)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqTransform"),
    .plotGenderMarkers.DESeqTransform)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("matrix"),
    .plotGenderMarkers)
