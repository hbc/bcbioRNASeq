#' @name params
#' @inherit basejump::params
#' @keywords internal
#'
#' @param normalized `character(1)`.
#'   Which normalization method to apply:
#'
#'     - "`tpm`": Transcripts per million (tximport).
#'     - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'     - "`rlog`": DESeq2 **log2** regularized log transformation.
#'     - "`vst`": DESeq2 **log2** variance stabilizing transformation.
#'
#'
#'
#' @param censorSamples `character`.
#'   Specify a subset of samples to censor.
#' @param color `ScaleDiscrete`.
#'   Desired ggplot2 color scale. Must supply discrete values. When set `NULL`,
#'   the default ggplot2 color palette will be used. If manual color definitions
#'   are desired, we recommend using [ggplot2::scale_color_manual()].
#'
#'   To set the discrete color palette globally, use:
#'
#'   ```
#'   options(acid.color.discrete = ggplot2::scale_color_viridis_d())
#'   ```
#' @param countsAxisLabel `character(1)`.
#'   Counts axis label.
#' @param fill `ggproto`/`ScaleDiscrete`.
#'   Desired ggplot2 fill scale. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#'
#'   To set the discrete fill palette globally, use:
#'
#'   ```
#'   options(acid.fill.discrete = ggplot2::scale_fill_viridis_d())
#'   ```
#' @param flip `logical(1)`.
#'   Flip x and y axes. Recommended for plots containing many samples.
#' @param geom `character(1)`.
#'   Plot type. Uses [`match.arg()`][base::match.arg] internally and defaults to
#'   the first argument in the `character` vector.
#' @param gffFile `character(1)`.
#'   GFF/GTF (General Feature Format) file. Generally, we recommend using a GTF
#'   (GFFv2) instead of a GFFv3 file if possible.
#' @param label `logical(1)`.
#'   Superimpose sample text labels on the plot.
#' @param legend `logical(1)`.
#'   Show plot legend.
#' @param limit `numeric(1)`.
#'   Threshold to denote on the plot, using a dashed line.
#' @param minCounts `integer(1)`.
#'   Minimum number of counts per gene in the count matrix.
#' @param perMillion `logical(1)`.
#'   Display as counts per million.
#' @param sampleMetadataFile `character(1)`.
#'   Sample metadata file path. CSV or TSV is preferred, but Excel worksheets
#'   are also supported. Check the documentation for conventions and required
#'   columns.
#' @param spikeNames `character`.
#'   Vector indicating which assay rows denote spike-in sequences (e.g. ERCCs).
#' @param title `character(1)`.
#'   Plot title.
#' @param trans `character(1)`.
#'   Name of the axis scale transformation to apply.
#'
#'   For more information:
#'
#'   ```
#'   help(topic = "scale_x_continuous", package = "ggplot2")
#'   ```
#' @param transgeneNames `character`.
#'   Vector indicating which assay rows denote transgenes (e.g. EGFP, TDTOMATO).
NULL
