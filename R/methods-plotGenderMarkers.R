#' Plot Sexually Dimorphic Gender Markers
#'
#' @name plotGenderMarkers
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotGene
#' @inheritParams general
#' @param organism Organism name. Typically can be left unset and should be
#'   detected automatically, unless a spike-in FASTA sequence is provided
#'   containing a gene identifier that is first alphabetically in the count
#'   matrix rownames.
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotGenderMarkers(bcb_small)
#'
#' # DESeqDataSet ====
#' plotGenderMarkers(dds_small)
#'
#' # DESeqTransform ====
#' plotGenderMarkers(rld_small)
NULL



# Constructors =================================================================
.plotGenderMarkers <- function(
    object,
    colData,
    interestingGroups = "sampleName",
    organism,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE
) {
    assert_is_matrix(object)
    assertFormalAnnotationCol(object, colData)
    assertFormalInterestingGroups(colData, interestingGroups)
    if (missing(organism)) {
        organism <- detectOrganism(object)
    }
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
    markers <- bcbioRNASeq::genderMarkers
    assert_is_subset(camel(organism), names(markers))
    markers <- markers[[camel(organism)]]

    # Ensembl identifiers
    gene2symbol <- markers %>%
        .[.[["include"]] == TRUE, , drop = FALSE] %>%
        mutate(
            geneName = paste(
                !!!syms(c("chromosome", "geneName")),
                sep = " : "
            )
        ) %>%
        .[, c("geneID", "geneName")] %>%
        as.data.frame() %>%
        set_rownames(.[["geneID"]])
    assertIsGene2symbol(gene2symbol)
    genes <- gene2symbol[["geneID"]] %>%
        .[. %in% rownames(object)]
    assert_is_non_empty(genes)

    .plotGene(
        object = object,
        genes = genes,
        gene2symbol = gene2symbol,
        colData = colData,
        interestingGroups = interestingGroups,
        countsAxisLabel = countsAxisLabel,
        medianLine = medianLine,
        color = color,
        return = "wide"
    ) +
        ggtitle(title)
}



# Methods ======================================================================
#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = c("rlog", "vst", "tpm"),
        color = scale_color_viridis(discrete = TRUE),
        title = TRUE
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        counts <- counts(object, normalized = normalized)
        # Ensure counts are log2 scale
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        .plotGenderMarkers(
            object = counts,
            interestingGroups = interestingGroups,
            organism = metadata(object)[["organism"]],
            colData = colData(object),
            countsAxisLabel = paste(normalized, "counts (log2)"),
            color = color,
            title = title
        )
    }
)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqDataSet"),
    function(
        object,
        interestingGroups = "sampleName",
        organism,
        color = scale_color_viridis(discrete = TRUE),
        title = TRUE
    ) {
        validObject(object)
        counts <- log2(counts(object, normalized = TRUE) + 1L)
        colData <- colData(object)
        # Drop the numeric sizeFactor column
        colData[["sizeFactor"]] <- NULL
        .plotGenderMarkers(
            object = counts,
            interestingGroups = interestingGroups,
            organism = organism,
            colData = colData,
            countsAxisLabel = "normalized counts (log2)",
            color = color,
            title = title
        )
    }
)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqTransform"),
    function(
        object,
        interestingGroups = "sampleName",
        organism,
        color = scale_color_viridis(discrete = TRUE),
        title = TRUE
    ) {
        # Passthrough: interestingGroups, color, title
        validObject(object)
        counts <- assay(object)

        colData <- colData(object)
        # Drop the numeric sizeFactor column
        colData[["sizeFactor"]] <- NULL
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        .plotGenderMarkers(
            object = counts,
            interestingGroups = interestingGroups,
            organism = organism,
            colData = colData,
            countsAxisLabel = paste(normalized, "counts (log2)"),
            color = color,
            title = title
        )
    }
)
