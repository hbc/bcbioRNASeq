#' Quality control metrics plots
#'
#' These functions are a collection of gene-level quality control plots.
#'
#' @rdname plot_metrics
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @param object Object.
#' @param pass_limit Threshold to plot pass color marker.
#' @param warn_limit Threshold to plot warning color marker.
#' @param interesting_groups *Optional*. Category to use to group samples (color
#'   and shape). If unset, this is automatically determined by the metadata set
#'   inside the [bcbioRNADataSet].
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#' @param filter_value Numeric value for filtering the counts matrix before
#'   plotting.
#'
#' @return [ggplot].
#'
#' @examples
#' data(bcb)
#' plot_total_reads(bcb)
#' plot_mapped_reads(bcb)
#' plot_mapping_rate(bcb)
#' plot_gene_detection_saturation(bcb)
#' plot_exonic_mapping_rate(bcb)
#' plot_intronic_mapping_rate(bcb)
#' plot_rrna_mapping_rate(bcb)
#' plot_counts_per_gene(bcb)
#' plot_count_density(bcb)



# Common functions ====
.assign_metrics_from_bcb <- function(bcb, envir = parent.frame()) {
    metrics(bcb) %>%
        assign("metrics", ., envir = envir)
    metadata(bcb)[["interesting_groups"]] %>%
        .[[1L]] %>%
        as.name %>%
        assign("interesting_group", ., envir = envir)
}



# [plot_total_reads()] ====
#' @rdname plot_metrics
#' @usage NULL
.plot_total_reads <- function(
    metrics,
    interesting_group,
    pass_limit = 20L,
    warn_limit = 10L,
    flip = TRUE) {
    p <- ggplot(metrics,
                aes_(x = ~description,
                     y = ~total_reads / 1e6L,
                     fill = as.name(interesting_group))) +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)",
             fill = "") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        geom_hline(alpha = 0.75,
                   color = pass_color,
                   size = 2L,
                   yintercept = pass_limit)
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot_metrics
#' @export
setMethod("plot_total_reads", "bcbioRNADataSet", function(
    object,
    pass_limit = 20L,
    warn_limit = 10L) {
    .assign_metrics_from_bcb(object)
    .plot_total_reads(metrics, interesting_group)
})







#' @rdname plot_metrics
#' @export
plot_mapped_reads <- function(
    bcb,
    pass_limit = 20L,
    warn_limit = 10L,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(
        metrics,
        aes_(x = ~description,
             y = ~mapped_reads / 1e6L,
             fill = interesting_groups)) +
        ggtitle("mapped reads") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        geom_hline(alpha = 0.75,
                   color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "mapped reads (million)",
             fill = "") +
        coord_flip()
}



#' @rdname plot_metrics
#' @export
plot_mapping_rate <- function(
    bcb,
    pass_limit = 90L,
    warn_limit = 70L,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(
        metrics,
        aes_(x = ~description,
             y = ~mapped_reads / total_reads * 100L,
             fill = interesting_groups)) +
        ggtitle("mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        geom_hline(alpha = 0.75,
                   color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}



#' @rdname plot_metrics
#' @export
plot_genes_detected <- function(
    bcb,
    pass_limit = 20000L,
    interesting_groups = NULL,
    filter_value = 0L) {
    metrics <- metrics(bcb)
    counts <- assay(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(metrics,
           aes_(x = ~description,
                y = colSums(counts > filter_value),
                fill = interesting_groups)) +
        ggtitle("genes detected") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "gene count",
             fill = "") +
        coord_flip()
}



#' @rdname plot_metrics
#' @export
plot_gene_detection_saturation <- function(
    bcb,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    counts <- counts(bcb)
    ggplot(
        metrics,
        aes_(x = ~mapped_reads / 1e6L,
             y = colSums(counts > 0L),
             color = interesting_groups)) +
        ggtitle("gene detection saturation") +
        geom_point(size = 3L) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(x = "mapped reads (million)",
             y = "gene count",
             color = "",
             shape = "")
}



#' @rdname plot_metrics
#' @export
plot_exonic_mapping_rate <- function(
    bcb,
    pass_limit = 60L,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(metrics,
           aes_(x = ~description,
                y = ~exonic_rate * 100L,
                fill = interesting_groups)) +
        ggtitle("exonic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "exonic mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}



#' @rdname plot_metrics
#' @export
plot_intronic_mapping_rate <- function(
    bcb,
    warn_limit = 20L,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(
        metrics,
        aes_(x = ~description,
             y = ~intronic_rate * 100L,
             fill = interesting_groups)) +
        ggtitle("intronic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "intronic mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}



#' @rdname plot_metrics
#' @export
plot_rrna_mapping_rate <- function(
    bcb,
    warn_limit = 10L,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(
        metrics,
        aes_(x = ~description,
             y = ~r_rna_rate * 100L,
             fill = interesting_groups)) +
        ggtitle("rRNA mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(alpha = 0.75,
                   color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "rRNA mapping rate (%)",
             fill = "") +
        coord_flip()
}



#' @rdname plot_metrics
#' @export
plot_counts_per_gene <- function(
    bcb,
    normalized = "tmm",
    interesting_groups = NULL) {
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    ggplot(
        melt_log10(bcb,
                   normalized = normalized,
                   interesting_groups = interesting_groups),
        aes_(x = ~description,
             y = ~counts,
             color = as.name(interesting_groups))) +
        ggtitle("counts per gene") +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "sample",
             y = "log10 counts per gene",
             color = interesting_groups) +
        coord_flip()
}


#' @rdname plot_metrics
#' @export
plot_count_density <- function(bcb, normalized = "tmm") {
    ggplot(
        melt_log10(bcb, normalized = normalized),
        aes_(x = ~counts,
             group = ~description)) +
        ggtitle("count density") +
        geom_density() +
        labs(x = "log10 counts per gene",
             y = "density")
}
