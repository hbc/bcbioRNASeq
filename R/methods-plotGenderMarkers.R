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
#' dir <- "http://bcbiornaseq.seq.cloud/f1000v1"
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
#' @importFrom ggplot2 aes_string element_text expand_limits geom_jitter ggplot
#'   labs theme
#' @importFrom magrittr set_colnames set_rownames
#' @importFrom tibble rownames_to_column
.plotGenderMarkers <- function(
    object,
    colData,
    interestingGroups = "sampleName",
    organism,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_matrix(object)
    assertFormalAnnotationCol(object, colData)
    assertFormalInterestingGroups(colData, interestingGroups)
    assert_is_a_string(organism)
    assert_is_a_string(countsAxisLabel)
    assert_is_a_bool(medianLine)
    assertIsColorScaleDiscreteOrNULL(color)

    # Title
    if (isTRUE(title)) {
        title <- "gender markers"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    colData <- uniteInterestingGroups(colData, interestingGroups)

    # Load the relevant internal gender markers data
    markers <- get("genderMarkers", envir = loadNamespace("bcbioRNASeq"))
    assert_is_subset(camel(organism), names(markers))
    markers <- markers[[camel(organism)]]

    # Ensembl identifiers
    gene2symbol <- markers %>%
        .[.[["include"]] == TRUE, , drop = FALSE] %>%
        mutate(
            symbol = paste(
                .data[["chromosome"]],
                .data[["symbol"]],
                sep = " : ")
        ) %>%
        .[, c("ensgene", "symbol")] %>%
        as.data.frame() %>%
        set_rownames(.[["ensgene"]])
    assertIsGene2symbol(gene2symbol)
    genes <- gene2symbol[["ensgene"]]
    assert_is_subset(genes, rownames(object))

    .plotGeneWide(
        object = object,
        genes = genes,
        gene2symbol = gene2symbol,
        colData = colData,
        interestingGroups = interestingGroups,
        countsAxisLabel = countsAxisLabel,
        medianLine = medianLine,
        color = color,
        title = title
    )
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
        colData = sampleMetadata(object),
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
        colData = sampleMetadata(object),
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
        colData = sampleMetadata(object),
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
