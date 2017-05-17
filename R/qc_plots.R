#' Standard RNA-seq quality control plots.
#'
#' @rdname qc_plots
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#' @author Victor Barrera
#'
#' @param run bcbio-nextgen run.
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   [DESeqDataSet] by running `counts(normalized = TRUE)`.
#'   Transcripts per million (TPM) are also acceptable.
#' @param raw_counts Raw counts matrix. Can be obtained from
#'   [DESeqDataSet] by running `counts(normalized = FALSE)`.
#' @param pass_limit Threshold to plot pass color marker.
#' @param warn_limit Threshold to plot warning color marker.
#'
#' @return ggplot2 plot object.



#' @rdname qc_plots
#' @export
plot_total_reads <- function(run, pass_limit = 20, warn_limit = 10) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        read_bcbio_metrics %>%
        ggplot(aes_(x = ~description,
                    y = ~total_reads / 1e6,
                    fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("total reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = warn_limit) +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "total reads (million)",
             fill = "") +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_mapped_reads <- function(run, pass_limit = 20, warn_limit = 10) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~description,
                 y = ~mapped_reads / 1e6,
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("mapped reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = warn_limit) +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "mapped reads (million)",
             fill = "") +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_mapping_rate <- function(run, pass_limit = 90, warn_limit = 70) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~description,
                 y = ~mapped_reads / total_reads * 100,
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = warn_limit) +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_genes_detected <- function(run, raw_counts, pass_limit = 20000) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~description,
                 y = colSums(raw_counts > 0),
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("genes detected") +
        geom_bar(stat = "identity") +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "gene count",
             fill = "") +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_gene_detection_saturation <- function(run, raw_counts) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~mapped_reads / 1e6,
                 y = ~colSums(raw_counts > 0),
                 color = as.name(run$intgroup[[1]]),
                 shape = as.name(run$intgroup[[1]]))) +
        ggtitle("gene detection saturation") +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(x = "mapped reads (million)",
             y = "gene count",
             color = "",
             shape = "")
}



#' @rdname qc_plots
#' @export
plot_exonic_mapping_rate <- function(run, pass_limit = 60) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~description,
                 # Multiple by 100 here for percentage
                 y = ~exonic_rate * 100,
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("exonic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = pass_limit) +
        labs(x = "sample",
             y = "exonic mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_intronic_mapping_rate <- function(run, warn_limit = 20) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~description,
                 # Multiple by 100 here for percentage
                 y = ~intronic_rate * 100,
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("intronic mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "intronic mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_rrna_mapping_rate <- function(run, warn_limit = 10) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    run$metrics %>%
        ggplot(
            aes_(x = ~description,
                 y = ~rrna_rate * 100,
                 fill = as.name(run$intgroup[[1]]))) +
        ggtitle("rRNA mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = warn_limit) +
        labs(x = "sample",
             y = "rRNA mapping rate (%)",
             fill = "") +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_counts_per_gene <- function(run, normalized_counts) {
    color <- run$intgroup[1]
    name <- deparse(substitute(normalized_counts))
    melted <- melt_log10(run, normalized_counts)
    data.frame(x = melted$description,
               y = melted$counts,
               color = melted[[color]]) %>%
        ggplot(
            aes_(x = ~x,
                 y = ~y,
                 color = ~color)
        ) +
        ggtitle(paste("counts per gene", name, sep = " : ")) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "sample",
             y = expression(log[10]~counts~per~gene),
             color = "") +
        coord_flip()
}



#' @rdname qc_plots
#' @export
plot_count_density <- function(
    run,
    normalized_counts) {
    name <- deparse(substitute(normalized_counts))
    melt_log10(run, normalized_counts) %>%
        ggplot(
            aes_(x = ~counts,
                 group = ~description)
        ) +
        ggtitle(paste("count density", name, sep = " : ")) +
        geom_density() +
        labs(x = expression(log[10]~counts~per~gene),
             y = "density")
}



#' @rdname qc_plots
#' @export
plot_gender_markers <- function(run, normalized_counts) {
    name <- deparse(substitute(normalized_counts))
    organism <- run$organism
    if (organism == "mmusculus") {
        gender_markers <- get("gender_markers_mmusculus",
                              envir = asNamespace("bcbioRnaseq"))
    } else {
        stop("Unsupported organism")
    }
    gender_markers <- gender_markers %>%
        .[.$include == TRUE, ] %>%
        .[order(.$chromosome, .$gene_symbol), ]
    # Ensembl identifiers
    identifier <- gender_markers$ensembl_gene %>% sort %>% unique
    # Scatterplot
    normalized_counts[identifier, ] %>%
        as.data.frame %>%
        rownames_to_column %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1) %>%
        set_names(c("ensembl_gene",
                    "description",
                    "counts")) %>%
        left_join(gender_markers, by = "ensembl_gene") %>%
        ggplot(
            aes_(x = ~gene_symbol,
                 y = ~counts,
                 color = ~description,
                 shape = ~chromosome)
        ) +
        ggtitle("gender markers") +
        geom_jitter(size = 4) +
        expand_limits(y = 0) +
        labs(x = "gene",
             y = name) +
        theme(legend.position = "none")
}
