#' Quality control plots
#'
#' @rdname qc_plots
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run \code{bcbio-nextgen} run
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   \code{DESeqDataSet} by running \code{counts(normalized = TRUE)}.
#'   Transcripts per million (TPM) are also acceptable.
#' @param raw_counts Raw counts matrix. Can be obtained from \code{DESeqDataSet}
#'   by running \code{counts(normalized = FALSE)}.



#' @rdname qc_plots
#' @description Plot total reads
#' @return Bar plot
#' @export
plot_total_reads <- function(run) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
        ggplot(aes_(x = ~description,
                    y = ~total_reads / 1e6,
                    fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("total reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = 10) +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = 20) +
        labs(x = "sample",
             y = "total reads (million)",
             fill = "") +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description Mapped reads plot
#' @return Bar plot
#' @export
plot_mapped_reads <- function(run) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
        ggplot(
            aes_(x = ~description,
                 y = ~mapped_reads / 1e6,
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("mapped reads") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = 10) +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = 20) +
        labs(x = "sample",
             y = "mapped reads (million)",
             fill = "") +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description Mapping rate plot
#' @return Bar plot
#' @export
plot_mapping_rate <- function(run) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
        ggplot(
            aes_(x = ~description,
                 y = ~mapped_reads / total_reads * 100,
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = 70) +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = 90) +
        labs(x = "sample",
             y = "mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description Genes detected plot
#' @return Bar plot
#' @export
plot_genes_detected <- function(run, raw_counts) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
        ggplot(
            aes_(x = ~description,
                 y = colSums(raw_counts > 0),
                 fill = as.name(run$intgroup[[1]]))
        ) +
        ggtitle("genes detected") +
        geom_bar(stat = "identity") +
        geom_hline(color = pass_color,
                   size = 2,
                   yintercept = 20000) +
        labs(x = "sample",
             y = "gene count",
             fill = "") +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description Genes detection saturation plot
#' @return Smooth plot
#' @export
plot_gene_detection_saturation <- function(run, raw_counts) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
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
    return(plot)
}



#' @rdname qc_plots
#' @description Exonic mapping rate plot
#' @return Bar plot
#' @export
plot_exonic_mapping_rate <- function(run) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
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
                   yintercept = 60) +
        labs(x = "sample",
             y = "exonic mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description Intronic mapping rate plot
#' @return Bar plot
#' @export
plot_intronic_mapping_rate <- function(run) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
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
                   yintercept = 20) +
        labs(x = "sample",
             y = "intronic mapping rate (%)",
             fill = "") +
        ylim(0, 100) +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description rRNA contamination mapping rate
#' @return Bar plot
#' @export
plot_rrna_mapping_rate <- function(run) {
    check_run(run)
    plot <- read_bcbio_metrics(run) %>%
        ggplot(
            aes_(x = ~description,
                 y = ~rrna_rate * 100,
                 fill = as.name(run$intgroup[[1]]))) +
        ggtitle("rRNA mapping rate") +
        geom_bar(stat = "identity") +
        geom_hline(color = warn_color,
                   size = 2,
                   yintercept = 10) +
        labs(x = "sample",
             y = "rRNA mapping rate (%)",
             fill = "") +
        coord_flip()
    return(plot)
}



#' @rdname qc_plots
#' @description Counts per gene plot
#' @return Boxplot
#' @export
plot_counts_per_gene <- function(run, normalized_counts) {
    check_run(run)
    color <- run$intgroup[1]
    name <- deparse(substitute(normalized_counts))
    melted <- melt_log10(run, normalized_counts)
    plot <- data.frame(x = melted$description,
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
    return(plot)
}



#' @rdname qc_plots
#' @description Plot count density
#' @return Density plot
#' @export
plot_count_density <- function(
    run,
    normalized_counts) {
    check_run(run)
    name <- deparse(substitute(normalized_counts))
    melted <- melt_log10(run, normalized_counts)
    plot <- ggplot(
        melted,
        aes_(x = ~counts,
             group = ~description)
    ) +
        ggtitle(paste("count density", name, sep = " : ")) +
        geom_density() +
        labs(x = expression(log[10]~counts~per~gene),
             y = "density")
    return(plot)
}



#' @rdname qc_plots
#' @description Plot sexually dimorphic gender marker genes
#' @return Scatter plot
#' @export
plot_gender_markers <- function(run, normalized_counts) {
    check_run(run)
    name <- deparse(substitute(normalized_counts))
    organism <- run$organism
    # Download the CSV file from seqcloud repo
    csv <- file.path("https://raw.githubusercontent.com",
                     "steinbaugh",
                     "seqcloud",
                     "master",
                     "workflows",
                     "illumina_rnaseq",
                     organism,
                     "gender_markers.csv") %>%
        read_csv(col_types = cols(),
                 na = c("", "#N/A")) %>%
        .[.$include == TRUE, ] %>%
        .[order(.$chromosome, .$gene_symbol), ]
    # Ensembl identifiers
    identifier <- csv[, "ensembl_gene"][[1]] %>% sort %>% unique
    # Scatterplot
    plot <- normalized_counts[identifier, ] %>%
        as.data.frame %>%
        rownames_to_column %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1) %>%
        set_names(c("ensembl_gene",
                    "description",
                    "counts")) %>%
        left_join(csv, by = "ensembl_gene") %>%
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
    return(plot)
}
