# [fix] Fast RNA-seq pipeline isn't outputting metrics. Does MultiQC 1.0 add
# support for lightweight/pseudoaligned counts?



#' Gene-level quality control plots
#'
#' @rdname qc_plots
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#' @author Victor Barrera
#'
#' @param run [bcbioRnaDataSet].
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   [DESeqDataSet] using [counts()] with `normalized = TRUE`.
#'   Transcripts per million (TPM) are also acceptable.
#' @param raw_counts Raw counts matrix. Can be obtained from
#'   [DESeqDataSet] using [counts()] with `normalized = FALSE`.
#' @param pass_limit Threshold to plot pass color marker.
#' @param warn_limit Threshold to plot warning color marker.
#' @param intgroup Category to use to group samples (color and shape). Default: `NULL`
#'
#' @return [ggplot].
#' @export
plot_total_reads <- function(run, pass_limit = 20, warn_limit = 10,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(total_reads = as.numeric(.data$total_reads)) %>%
        ggplot(aes_(x = ~description,
                    y = ~total_reads / 1e6,
                    fill = intgroup)) +
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
plot_mapped_reads <- function(run, pass_limit = 20, warn_limit = 10,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(mapped_reads = as.numeric(.data$mapped_reads)) %>%
        ggplot(aes_(x = ~description,
                    y = ~mapped_reads / 1e6,
                    fill = intgroup)) +
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
plot_mapping_rate <- function(run, pass_limit = 90, warn_limit = 70,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(mapped_reads = as.numeric(.data$mapped_reads),
               total_reads = as.numeric(.data$total_reads)) %>%
        ggplot(aes_(x = ~description,
                    y = ~mapped_reads / total_reads * 100,
                    fill = intgroup)) +
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
plot_genes_detected <- function(run, raw_counts, pass_limit = 20000,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        ggplot(aes_(x = ~description,
                    y = colSums(raw_counts > 0),
                    fill = intgroup)) +
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
plot_gene_detection_saturation <- function(run, raw_counts,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(mapped_reads = as.numeric(.data$mapped_reads)) %>%
        ggplot(aes_(x = ~mapped_reads / 1e6,
                    y = ~colSums(raw_counts > 0),
                    color = intgroup,
                    shape = intgroup)) +
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
plot_exonic_mapping_rate <- function(run, pass_limit = 60,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(exonic_rate = as.numeric(.data$exonic_rate)) %>%
        ggplot(aes_(x = ~description,
                    y = ~exonic_rate * 100,
                    fill = intgroup)) +
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
plot_intronic_mapping_rate <- function(run, warn_limit = 20,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(intronic_rate = as.numeric(.data$intronic_rate)) %>%
        ggplot(aes_(x = ~description,
                    y = ~intronic_rate * 100,
                    fill = intgroup)) +
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
plot_rrna_mapping_rate <- function(run, warn_limit = 10,intgroup=NULL) {
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
    }
    run$metrics %>%
        mutate(r_rna_rate = as.numeric(.data$r_rna_rate) * 100) %>%
        ggplot(aes_(x = ~description,
                    y = ~r_rna_rate,
                    fill = intgroup)) +
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
plot_counts_per_gene <- function(run, normalized_counts,intgroup=NULL) {
    name <- deparse(substitute(normalized_counts))
    melted <- melt_log10(run, normalized_counts)

    if(!is.null(intgroup)){
        intgroup=as.name(intgroup)
        melted$color=melted[[intgroup]]
    }else{
        melted$color="ALL"
    }

    data.frame(x = melted$description,
               y = melted$counts,
               color = melted$color) %>%
        ggplot(aes_(x = ~x,
                    y = ~y,
                    color = ~color)) +
        ggtitle(paste("counts per gene", name, sep = label_sep)) +
        geom_boxplot(outlier.shape = NA) +
        # optional way to make log10 subscript:
        # expression(log[10]~counts~per~gene)
        labs(x = "sample",
             y = "log10 counts per gene",
             color = intgroup) +
        coord_flip()
}


#' @rdname qc_plots
#' @export
plot_count_density <- function(
    run,
    normalized_counts) {
    name <- deparse(substitute(normalized_counts))
    melt_log10(run, normalized_counts) %>%
        ggplot(aes_(x = ~counts,
                    group = ~description)) +
        ggtitle(paste("count density", name, sep = label_sep)) +
        geom_density() +
        labs(x = "log10 counts per gene",
             y = "density")
}



#' @rdname qc_plots
#' @export
plot_gender_markers <- function(run, normalized_counts) {
    name <- deparse(substitute(normalized_counts))
    organism <- run$organism

    # Load the relevant internal gender markers data
    if (organism == "mmusculus") {
        gender_markers <- get("gender_markers_mmusculus",
                              envir = asNamespace("bcbioRnaseq"))
    } else {
        stop("Unsupported organism")
    }

    # Prepare the source gender markers data frame
    gender_markers <- gender_markers %>%
        filter(.data$include == TRUE) %>%
        arrange(!!!syms(c("chromosome", "gene_symbol")))

    # Ensembl identifiers
    identifier <- gender_markers$ensembl_gene %>% sort %>% unique

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



#' Find correlation between PCs and covariates
#'
#' @author Lorena Pantano
#'
#' @param run [bcbioRnaDataSet].
#' @param dt [DESeqTransform]. [rlog()]-transformed counts are recommended.
#' @param use character vector. List of columns to use in degCovariates
#' @param ... Passthrough parmeters to [DEGreport::degCovariates()].
#'
#' @export
plot_pca_covariates <- function(run, dt, use=NULL, ...) {
    # [fix] Not working for some consults?
    if (is.null(run$metrics)) {
        return(NULL)
    }
    if (is.null(use)){
        use <- colnames(run$metrics)
    }else{
        use <- intersect(use, colnames(run$metrics))
    }
    if (length(use) == 0)
        stop("Not columns matched between use and metadata")

    keep_metrics <- lapply(use, function(a) {
        if (length(unique(run$metrics[, a])) > 1) { a }
    }) %>% unlist %>% .[!is.null(.)]
    metrics <- run$metrics %>%
        as.data.frame %>%
        set_rownames(.$description) %>%
        .[, setdiff(keep_metrics, c("description", "file_name"))]
    dt %>%
        assay %>%
        degCovariates(metadata = metrics, ...) %>%
        .$plot +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
}
