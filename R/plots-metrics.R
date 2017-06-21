#' Quality control metrics plots
#'
#' These functions are a collection of gene-level quality control plots.
#'
#' @rdname plots-metrics
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @param bcb [bcbioRnaDataSet].
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#' @param pass_limit Threshold to plot pass color marker.
#' @param warn_limit Threshold to plot warning color marker.
#' @param interesting_groups (*Optional*). Category to use to group samples
#'   (color and shape). If unset, this is automatically determined by
#'   the metadata set inside the [bcbioRnaDataSet].
#'
#' @return [ggplot].
#' @export
plot_total_reads <- function(
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
             y = ~total_reads / 1e6,
             fill = interesting_groups)) +
        ggtitle("total reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        geom_hline(color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "total reads (million)",
             fill = "") +
        coord_flip()
}



#' @rdname plots-metrics
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
    ggplot(
        metrics,
        aes_(x = ~description,
             y = ~mapped_reads / 1e6,
             fill = interesting_groups)) +
        ggtitle("mapped reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        geom_hline(color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "mapped reads (million)",
             fill = "") +
        coord_flip()
}



#' @rdname plots-metrics
#' @export
plot_mapping_rate <- function(
    bcb,
    pass_limit = 90,
    warn_limit = 70,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    ggplot(
        metrics,
        aes_(x = ~description,
             y = ~mapped_reads / total_reads * 100L,
             fill = interesting_groups)) +
        ggtitle("mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        geom_hline(color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}



#' @rdname plots-metrics
#' @export
plot_genes_detected <- function(
    bcb,
    pass_limit = 20000,
    interesting_groups = NULL) {
    metrics <- metrics(bcb)
    counts <- assay(bcb)
    if (is.null(metrics)) return(NULL)
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }
    interesting_groups <- as.name(interesting_groups)
    ggplot(metrics,
           aes_(x = ~description,
                y = colSums(counts > 0L),
                fill = interesting_groups)) +
        ggtitle("genes detected") +
        geom_bar(stat = "identity") +
        geom_hline(color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "gene count",
             fill = "") +
        coord_flip()
}



#' @rdname plots-metrics
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
        aes_(x = ~mapped_reads / 1e6,
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



#' @rdname plots-metrics
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
        geom_hline(color = pass_color,
                   size = 2L,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "exonic mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}



#' @rdname plots-metrics
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
        geom_hline(color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "intronic mapping rate (%)",
             fill = "") +
        ylim(0L, 100L) +
        coord_flip()
}



#' @rdname plots-metrics
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
        geom_hline(color = warn_color,
                   size = 2L,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "rRNA mapping rate (%)",
             fill = "") +
        coord_flip()
}



#' @rdname plots-metrics
#' @export
plot_counts_per_gene <- function(
    bcb,
    counts = "tmm",
    interesting_groups = NULL) {
    if (counts == "tmm") {
        title <- "counts per gene (tmm normalized)"
    } else {
        title <- "counts per gene"
    }

    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(bcb)[["interesting_groups"]]
    }

    ggplot(
        melt_log10(bcb, counts),
        aes_(x = ~description,
             y = ~counts,
             color = as.name(interesting_groups))) +
        ggtitle(title) +
        geom_boxplot(outlier.shape = NA) +
        # optional way to make log10 subscript:
        # expression(log[10]~counts~per~gene)
        labs(x = "sample",
             y = "log10 counts per gene",
             color = interesting_groups) +
        coord_flip()
}


#' @rdname plots-metrics
#' @export
plot_count_density <- function(bcb, counts = "tmm") {
    if (counts == "tmm") {
        title <- "counts per gene (tmm-normalized)"
    } else {
        title <- "counts per gene"
    }
    ggplot(
        melt_log10(bcb, counts),
        aes_(x = ~counts,
             group = ~description)) +
        ggtitle(title) +
        geom_density() +
        labs(x = "log10 counts per gene",
             y = "density")
}
