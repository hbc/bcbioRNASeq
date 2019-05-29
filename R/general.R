#' General arguments
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
#' @param alpha `scalar numeric`. Adjusted P value ("alpha") cutoff.
#' @param color `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 color
#'   scale. Must supply discrete values. When set to `NULL`, the default ggplot2
#'   color palette will be used. If manual color definitions are desired, we
#'   recommend using [ggplot2::scale_color_manual()].
#'   To set the discrete color palette globally, use
#'   `options(bcbio.discrete.color = scale_color_viridis_d())`.
#' @param counts `matrix`. Normalized counts.
#' @param dir `string`. Local directory path.
#' @param direction `string`. Plot "`both`", "`up`", or "`down`" directions.
#' @param fill `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 fill scale.
#'   Must supply discrete values. When set to `NULL`, the default ggplot2 color
#'   palette will be used. If manual color definitions are desired, we recommend
#'   using [ggplot2::scale_fill_manual()].
#'   To set the discrete fill palette globally, use
#'   `options(bcbio.discrete.fill = scale_fill_viridis_d())`.
#' @param flip `boolean`. Flip x and y axes. Recommended for quality control
#'   plots containing many samples.
#' @param gene2symbol `data.frame`. Gene-to-symbol mappings. Must contain the
#'   columns `geneID` and `geneName`.
#' @param genes `character`. Genes to include. These must match the rownames of
#'   the object. It is best practice to use the stable gene identifiers from
#'   Ensembl (e.g. "ENSG00000000003") and not the gene symbols.
#' @param headerLevel `scalar integer` (`1`-`7`). R Markdown header level.
#' @param interestingGroups `character`. Groups of interest that define the
#'   samples. If left unset, defaults to `sampleName`.
#' @param label `boolean`. Superimpose sample text labels on the plot.
#' @param legend `boolean`. Show plot legend.
#' @param lfcThreshold `scalar numeric`. Log fold change ratio (base 2) cutoff
#'   threshold.
#' @param limit `scalar numeric`. Threshold to denote on the plot, using a
#'   dashed line.
#' @param minCounts `scalar integer`. Minimum number of counts per gene in the
#'   counts matrix.
#' @param normalized `string`. Which normalization method to apply:
#'   - "`tpm`": Transcripts per million (tximport).
#'   - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'   - "`rlog`": DESeq2 **log2** regularized log transformation.
#'   - "`vst`": DESeq2 **log2** variance stabilizing transformation.
#' @param ntop `scalar integer`. Number of top genes to label.
#' @param passLimit `scalar numeric`. Threshold to plot pass color marker.
#' @param pointColor `string`. Default point color for the plot.
#' @param results `DESeqResults`.
#' @param return `string`. Object class to return. Uses [match.arg()] internally
#'   and picks the first item in the vector by default.
#' @param samples `character` or `NULL`. Samples to include.
#' @param sigPointColor `character`. Color names for labeling upregulated and
#'   downregulated genes. Also supports a character string for labeling DEGs
#'   with the same color, regardless of direction.
#' @param subtitle `string` or `NULL`. Plot subtitle.
#' @param title `string` or `NULL`. Plot title.
#' @param warnLimit `scalar numeric`. Threshold to plot warning color marker.
#' @param withDimnames `boolean`. Whether dimnames should be applied to
#'   extracted assay elements.
#'
#' @return No value.
NULL
