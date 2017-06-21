#' Gene-level quality control plots
#'
#' @rdname qc_plots
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#' @author Victor Barrera
#'
#' @param bcb [bcbioRnaDataSet].
#' @param normalized Apply normalization to counts.
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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



#' @rdname qc_plots
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


#' @rdname qc_plots
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



#' @rdname qc_plots
#' @export
plot_gender_markers <- function(bcb, normalized_counts = NULL) {
    # Organism-specific dimorphic markers ----
    organism <- metadata(bcb)[["organism"]]

    # Load the relevant internal gender markers data
    if (organism == "mmusculus") {
        gender_markers <- get("gender_markers_mmusculus",
                              envir = asNamespace("bcbioRnaseq"))
    } else {
        stop("Unsupported organism")
    }

    # Prepare the source gender markers data frame
    gender_markers <- gender_markers %>%
        filter(.data[["include"]] == TRUE) %>%
        arrange(!!!syms(c("chromosome", "gene_symbol")))

    # Ensembl identifiers
    identifier <- gender_markers[["ensembl_gene"]] %>%
        sort %>% unique


    # Normalized counts ----
    if (is.null(normalized_counts)) {
        normalized_counts <- tpm(bcb)
        name <- "transcripts per million (tpm)"
    } else {
        name <- deparse(substitute(normalized_counts))
    }


    # ggplot ----
    normalized_counts[identifier, ] %>%
        as.data.frame %>%
        rownames_to_column %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1) %>%
        set_names(c("ensgene",
                    "description",
                    "counts")) %>%
        left_join(gender_markers, by = "ensgene") %>%
        ggplot(
            aes_(x = ~gene_symbol,
                 y = ~counts,
                 color = ~description,
                 shape = ~chromosome)) +
        ggtitle("gender markers") +
        geom_jitter(size = 4) +
        expand_limits(y = 0) +
        labs(x = "gene",
             y = name)
}



#' Find correlation between PCs and covariates
#'
#' @author Lorena Pantano
#'
#' @param bcb [bcbioRnaDataSet].
#' @param dt [DESeqTransform]. [rlog()]-transformed counts are recommended.
#' @param use Character vector. List of columns to use in degCovariates.
#' @param ... Passthrough parmeters to [DEGreport::degCovariates()].
#'
#' @export
plot_pca_covariates <- function(bcb, dt, use = NULL, ...) {
    metrics <- metrics(bcb)
    if (is.null(metrics)) {
        return(NULL)
    }

    if (is.null(use)) {
        use <- colnames(metrics)
    } else {
        use <- intersect(use, colnames(metrics))
    }

    if (length(use) == 0L) {
        stop("Not columns matched between use and metadata")
    }

    keep_metrics <- lapply(use, function(a) {
        if (length(unique(metrics[, a])) > 1L) a
    }) %>% unlist %>% .[!is.null(.)]

    metrics <- metrics %>%
        as.data.frame %>%
        set_rownames(.$description) %>%
        .[, setdiff(keep_metrics, c("description", "file_name"))]

    dt %>%
        assay %>%
        degCovariates(metadata = metrics, ...) %>%
        .[["plot"]] +
        theme(axis.text.x = element_text(angle = 60L, hjust = 1L))
}
