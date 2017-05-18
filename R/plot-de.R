#' Differential expression plots
#'
#' @rdname de_plots
#' @author Michael Steinbaugh
#'
#' @seealso [DESeq2::plotMA()].



#' @rdname de_plots
#' @description Wrapper for [DESeq2::plotMA()] that generates a title
#'   automatically.
#'
#' @param res \linkS4class{DESeqResults}.
#' @param ylim Y-axis maximum (single integer).
#'
#' @return MA plot.
#' @export
plot_ma <- function(res, ylim = 2) {
    check_res(res)

    name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)

    plotMA(
        res,
        main = paste(name, contrast_name, sep = " : "),
        ylim = c(-ylim, ylim))
}



# Modified version of `CHBUtils::volcano_density_plot()`
# https://github.com/hbc/CHBUtils/blob/master/R/volcanoPlot.R

#' @rdname de_plots
#' @description Volcano plot.
#'
#' @param run bcbio-nextgen run.
#' @param lfc Log fold change ratio (base 2) cutoff for coloring.
#' @param text_labels Number of text labels to plot.
#'
#' @param direction Plot `up`, `down`, or `both` directions.
#' @param title Title for the figure.
#' @param shade_color Shading color for bounding box.
#' @param shade_alpha Shading transparency alpha.
#' @param point_color Point color.
#' @param point_alpha Point transparency alpha.
#' @param point_outline_color Point outline color.
#'
#' @return Volcano plot.
#' @export
plot_volcano <- function(
    run,
    res,
    lfc = 1,
    text_labels = 30,
    direction = "both",
    title = NULL,
    shade_color = "green",
    shade_alpha = 0.25,
    point_color = "gray",
    point_alpha = 0.75,
    point_outline_color = "darkgray") {
    check_run(run)
    check_res(res)
    import_tidy_verbs()

    if (!any(direction %in% c("both", "down", "up")) |
        length(direction) > 1) {
        stop("direction must be both, up, or down")
    }

    alpha <- res@metadata$alpha
    contrast_name <- res_contrast_name(res)

    # Generate automatic title, if necessary
    if (is.null(title)) {
        title <- paste(deparse(substitute(res)),
                       contrast_name,
                       sep = " : ")
    }

    stats <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        snake %>%
        select(!!!syms(c("ensembl_gene_id", "log2_fold_change", "padj"))) %>%
        # Filter zero counts for quicker plotting
        filter(!is.na(.data$log2_fold_change)) %>%
        # Arrange by P value
        arrange(!!sym("padj")) %>%
        mutate(neg_log_padj = -log10(.data$padj + 1e-10))

    # Automatically label the top genes
    annotations <- gene_level_annotations(run) %>%
        select(!!!syms(c("ensembl_gene_id", "external_gene_name")))
    volcano_text <- stats %>%
        .[1:text_labels, ] %>%
        left_join(annotations, by = "ensembl_gene_id")

    # Get range of LFC and P values to set up plot borders
    range_lfc <-
        c(floor(min(na.omit(stats$log2_fold_change))),
          ceiling(max(na.omit(stats$log2_fold_change))))
    range_neg_log_padj <-
        c(floor(min(na.omit(stats$neg_log_padj))),
          ceiling(max(na.omit(stats$neg_log_padj))))


    # LFC density histogram ====
    lfc_density <- density(na.omit(stats$log2_fold_change))
    lfc_density_df <- data.frame(x = lfc_density$x,
                                 y = lfc_density$y)
    lfc_hist <- stats %>%
        ggplot(aes_(x = ~log2_fold_change)) +
        geom_density() +
        scale_x_continuous(limits = range_lfc) +
        labs(title = "log2 fold change density plot",
             x = "log2 fold change")
    if (direction == "both" | direction == "up") {
        lfc_hist <- lfc_hist +
            geom_ribbon(
                data = filter(lfc_density_df,
                              .data$x > UQ(lfc)),
                aes_(x = ~x, ymax = ~y),
                ymin = 0,
                fill = shade_color,
                alpha = shade_alpha)
    }
    if (direction == "both" | direction == "down") {
        lfc_hist <- lfc_hist +
            geom_ribbon(
                data = filter(lfc_density_df,
                              .data$x < -UQ(lfc)),
                aes_(x = ~x, ymax = ~y),
                ymin = 0,
                fill = shade_color,
                alpha = shade_alpha)
    }
    show(lfc_hist)


    # Density plot of adjusted P values ====
    padj_density <- density(na.omit(stats$neg_log_padj))
    padj_density_df <- data.frame(x = padj_density$x,
                                  y = padj_density$y)
    padj_hist <- stats %>%
        ggplot(aes_(x = ~-log10(padj + 1e-10))) +
        geom_density() +
        # filter(lfc_density_df, .data$x < -lfc),
        geom_ribbon(data = filter(padj_density_df,
                                  .data$x > -log10(UQ(alpha) + 1e-10)),
                    aes_(x = ~x, ymax = ~y),
                    ymin = 0,
                    fill = shade_color,
                    alpha = shade_alpha) +
        coord_flip() +
        labs(title = "p value density plot",
             y = "density")
    show(padj_hist)


    # Volcano plot ====
    volcano <- stats %>%
        ggplot(aes_(x = ~log2_fold_change,
                    y = ~-log10(padj + 1e-10))) +
        labs(title = "volcano plot",
             x = "log2 fold change") +
        geom_point(
            alpha = point_alpha,
            color = point_outline_color,
            fill = point_color,
            pch = 21) +
        theme(legend.position = "none") +
        scale_x_continuous(limits = range_lfc) +
        geom_text_repel(
            data = volcano_text,
            aes_(x = ~log2_fold_change,
                 y = ~neg_log_padj,
                 label = ~external_gene_name),
            size = 3)
    if (direction == "both" | direction == "up") {
        volcano_poly_up <- with(stats, data.frame(
            x = as.numeric(c(
                lfc,
                lfc,
                max(range_lfc),
                max(range_lfc))),
            y = as.numeric(c(
                -log10(alpha + 1e-10),
                max(range_neg_log_padj),
                max(range_neg_log_padj),
                -log10(alpha + 1e-10)))))
        volcano <- volcano +
            geom_polygon(
                data = volcano_poly_up,
                aes_(x = ~x, y = ~y),
                fill = shade_color,
                alpha = shade_alpha)
    }
    if (direction == "both" | direction == "down") {
        volcano_poly_down <- with(stats, data.frame(
            x = as.numeric(c(
                -lfc,
                -lfc,
                min(range_lfc),
                min(range_lfc))),
            y = as.numeric(c(
                -log10(alpha + 1e-10),
                max(range_neg_log_padj),
                max(range_neg_log_padj),
                -log10(alpha + 1e-10)))))
        volcano <- volcano +
            geom_polygon(
                data = volcano_poly_down,
                aes_(x = ~x, y = ~y),
                fill = shade_color,
                alpha = shade_alpha)
    }
    show(volcano)
}



#' Differentially expressed gene heatmap
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param res \linkS4class{DESeqResults}.
#' @param dt \linkS4class{DESeqTransform}. We recommend passing in [rlog()]
#'   counts in as the transform object by default.
#' @param lfc Log2 fold change ratio cutoff.
#'
#' @return Gene heatmap.
#' @export
plot_deg_heatmap <- function(
    run,
    res,
    dt,
    lfc = 0) {
    check_run(run)
    check_res(res)
    check_dt(dt)
    import_tidy_verbs()

    alpha <- res@metadata$alpha
    metadata <- intgroup_as_factor(run)

    # Output results to data frame and subset by alpha and lfc cutoffs
    res_df <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        filter(.data$padj < alpha) %>%
        .[which(.$log2FoldChange > lfc | .$log2FoldChange < -lfc), ]

    # rlog transform and subset the counts
    counts <- dt %>% assay %>% .[res_df$ensembl_gene_id, ]

    if (nrow(counts) <= 100) {
        # Change rownames to readable external gene names
        show_rownames <- TRUE
        gene_names <- gene_level_annotations(run) %>%
            select(!!!syms(c("ensembl_gene_id", "external_gene_name"))) %>%
            as.data.frame %>%
            set_rownames(.$ensembl_gene_id)
        rownames(counts) <- gene_names[rownames(counts), "external_gene_name"]
    } else {
        show_rownames <- FALSE
    }

    res_name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)
    dt_name <- deparse(substitute(dt))

    pheatmap(counts,
             annotation = metadata,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D2",
             main = paste(res_name, contrast_name, dt_name, sep = " : "),
             scale = "row",
             show_rownames = show_rownames)
}
