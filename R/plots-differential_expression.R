#' Differential expression plots
#'
#' @rdname de_plots
#' @author Michael Steinbaugh



#' @rdname de_plots
#' @description Wrapper for [DESeq2::plotMA()] that generates a title
#'   automatically.
#'
#' @param res [DESeqResults].
#' @param ylim Y-axis maximum (single integer).
#'
#' @return MA plot.
#' @export
plot_ma <- function(res, ylim = 2) {
    .check_res(res)
    name <- deparse(substitute(res))
    contrast_name <- .res_contrast_name(res)
    plotMA(res,
           main = paste(name, contrast_name, sep = label_sep),
           ylim = c(-ylim, ylim))
}



#' @rdname de_plots
#' @description Volcano plot.
#'
#' @param bcb [bcbioRnaDataSet].
#' @param lfc Log fold change ratio (base 2) cutoff for coloring.
#' @param text_labels Number of text labels to plot.
#' @param merge merge all plots into one.
#' @param direction Plot `up`, `down`, or `both` directions.
#' @param title Title for the figure.
#' @param shade_color Shading color for bounding box.
#' @param shade_alpha Shading transparency alpha.
#' @param point_color Point color.
#' @param point_alpha Point transparency alpha.
#' @param point_outline_color Point outline color.
#'
#' @author Michael Steinbaugh and Lorena Pantano (based on John Hutchinson's
#'   work)
#' @return Volcano plot.
#' @export
#'
#' @seealso This function is a modified version of
#'   `CHBUtils::volcano_density_plot()`.
plot_volcano <- function(
    bcb,
    res,
    lfc = 1L,
    text_labels = 30L,
    merge = TRUE,
    direction = "both",
    title = NULL,
    shade_color = "orange",
    shade_alpha = 0.25,
    point_color = "gray",
    point_alpha = 0.75,
    point_outline_color = "darkgray") {
    .check_res(res)
    # [TODO] Add support for option of plotting unadjusted P values without
    # `+ 1e-10` transformation.

    if (!any(direction %in% c("both", "down", "up")) |
        length(direction) > 1) {
        stop("direction must be both, up, or down")
    }

    alpha <- metadata(res)$alpha
    contrast_name <- .res_contrast_name(res)

    # Generate automatic title, if necessary
    if (is.null(title)) {
        title <- paste(deparse(substitute(res)),
                       contrast_name,
                       sep = label_sep)
    }

    stats <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        snake %>%
        tidy_select(!!!syms(c("ensembl_gene_id", "log2_fold_change", "padj"))) %>%
        # Filter zero counts for quicker plotting
        filter(!is.na(.data$log2_fold_change)) %>%
        # Arrange by P value
        arrange(!!sym("padj")) %>%
        # Convert adjusted P value to -log10
        mutate(neg_log_padj = -log10(.data$padj + 1e-10))

    # Automatically label the top genes
    # [TODO] Make annotable interconversion to Ensembl names a function...
    annotations <- rowData(bcb)[, c("ensgene", "symbol")] %>%
        as.data.frame %>%
        set_colnames(c("ensembl_gene_id", "external_gene_name"))

    volcano_text <- stats[1:text_labels, ] %>%
        left_join(annotations, by = "ensembl_gene_id")

    # Get range of LFC and P values to set up plot borders
    range_lfc <-
        c(floor(min(na.omit(stats$log2_fold_change))),
          ceiling(max(na.omit(stats$log2_fold_change))))
    range_neg_log_padj <-
        c(floor(min(na.omit(stats$neg_log_padj))),
          ceiling(max(na.omit(stats$neg_log_padj))))


    # LFC density histogram ----
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
                # [fix] check UQ syntax
                data = filter(lfc_density_df, .data$x > UQ(lfc)),
                aes_(x = ~x, ymax = ~y),
                ymin = 0,
                fill = shade_color,
                alpha = shade_alpha)
    }
    if (direction == "both" | direction == "down") {
        lfc_hist <- lfc_hist +
            geom_ribbon(
                # [fix] check UQ syntax
                data = filter(lfc_density_df, .data$x < -UQ(lfc)),
                aes_(x = ~x, ymax = ~y),
                ymin = 0,
                fill = shade_color,
                alpha = shade_alpha)
    }
    # show(lfc_hist)


    # Density plot of adjusted P values ----
    padj_density <- density(na.omit(stats$neg_log_padj))
    padj_density_df <- data.frame(x = padj_density$x,
                                  y = padj_density$y)
    padj_hist <- stats %>%
        ggplot(aes_(x = ~-log10(padj + 1e-10))) +
        geom_density() +
        # [fix] check that `UQ` is best approach here
        geom_ribbon(data = filter(padj_density_df,
                                  .data$x > -log10(UQ(alpha) + 1e-10)),
                    aes_(x = ~x, ymax = ~y),
                    ymin = 0,
                    fill = shade_color,
                    alpha = shade_alpha) +
        coord_flip() +
        labs(title = "p-value density plot",
             y = "density")
    # show(padj_hist)


    # Volcano plot ----
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


    # Grid layout ----
    if (merge) {
        ggdraw() +
            draw_plot(
                lfc_hist +
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank()),
                x = 0L, y = 0.7, width = 1L, height = 0.3) +
            draw_plot(
                padj_hist +
                    theme(axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank()),
                x = 0.7, y = 0L, width = 0.3, height = 0.7) +
            draw_plot(
                volcano, x = 0L, y = 0L, width = 0.7, height = 0.7)
    } else {
        show(lfc_hist)
        show(padj_hist)
        show(volcano)
    }
}


#' Wrap degPatter function around bcbio object
#'
#' It runs [DEGreport::degPatterns] function from a bcbio object
#'
#' @param bcb [bcbioRnaseq] object
#' @param res table with padj as column and gene names as row.names
#' @param fdr float cutoff to consider genes significant
#' @param n integer maximum number of genes to consider
#' @param ... extra parameters for degPatterns function
#' @author Lorena Pantano
#' @export
plot_pattern <- function(bcb, res, fdr = 0.1, n = NULL, ...){
    res <- as.data.frame(res)
    .order <- res[order(res$padj), ]
    sign <- row.names(.order)[!is.na(.order$padj)]
    if (!is.null(n))
        sign <- sign[1:n]
    degPatterns(bcbio(bcb, "rld")[sign, ], colData(bcb), ...)
}

