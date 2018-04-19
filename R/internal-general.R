#' General Arguments
#'
#' @name general
#' @keywords internal
#'
#' @param object Object.
#' @param i An integer or numeric scalar.
#' @param value Value to assign.
#' @param x Object.
#' @param ... Additional arguments.
#'
#' @param alpha Adjusted P value ("alpha") cutoff.
#' @param colData `data.frame` describing the columns of the object. Must
#'   contain `rownames` identical to the `colnames` of the object.
#' @param color Desired ggplot color scale. Must supply discrete values. When
#'   set to `NULL`, the default ggplot2 color palette will be used. If manual
#'   color definitions are desired, we recommend using
#'   [ggplot2::scale_color_manual()].
#' @param counts Object containing a normalized counts matrix.
#' @param dir Local directory path.
#' @param fill Desired ggplot fill scale. Must supply discrete values. When set
#'   to `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#' @param flip Flip x and y axes.
#' @param gene2symbol `data.frame` containing gene-to-symbol mappings. Must
#'   contain the columns `geneID` and `geneName`.
#' @param genes Character vector of genes to include. These must match the
#'   rownames of the object. It is best practice to use the stable gene
#'   identifiers from Ensembl (e.g. "ENSG00000000003") and not the gene symbols.
#' @param headerLevel R Markdown header level.
#' @param interestingGroups Character vector denoting groups of interest that
#'   define the samples. If left unset, defaults to `sampleName`.
#' @param lfcThreshold Log fold change ratio (base 2) cutoff threshold.
#' @param minCounts Minimum number of counts per gene in the counts matrix.
#' @param normalized Character indicating which normalization method to apply:
#'   - "`tpm`": Transcripts per million (tximport).
#'   - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'   - "`rlog`": DESeq2 **log2** regularized log transformation.
#'   - "`vst`": DESeq2 **log2** variance stabilizing transformation.
#' @param ntop Number of top genes to label.
#' @param passLimit Threshold to plot pass color marker.
#' @param pointColor Default point color for the plot.
#' @param results `DESeqResults` object.
#' @param return Object class to return. Uses [match.arg()] internally and picks
#'   the first item in the vector by default.
#' @param samples Character vector of samples to include.
#' @param sigPointColor `character` vector containing color names for labeling
#'   upregulated and downregulated genes. Also supports a character string for
#'   labeling DEGs with the same color, regardless of direction.
#' @param subtitle Subtitle of plot.
#' @param title Title of plot.
#' @param warnLimit Threshold to plot warning color marker.
#' @param withDimnames A `logical`, indicating whether dimnames should be
#'   applied to extracted assay elements.
#'
#' @return No value.
NULL
