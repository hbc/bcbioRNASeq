#' Plot Volcano
#'
#' @rdname plot_volcano
#' @author John Hutchinson, Michael Steinbaugh, Lorena Pantano
#' @family Differential Expression Plots
#'
#' @param alpha Alpha level cutoff used for coloring.
#' @param padj Use P values adjusted for multiple comparisions.
#' @param lfc Log fold change ratio (base 2) cutoff for coloring.
#' @param genes Character vector of gene symbols to label.
#' @param ntop Number of top genes to label.
#' @param direction Plot `up`, `down`, or `both` (**default**) directions.
#' @param shade_color Shading color for bounding box.
#' @param shade_alpha Shading transparency alpha.
#' @param point_color Point color.
#' @param point_alpha Point transparency alpha.
#' @param point_outline_color Point outline color.
#' @param grid Arrange plots into grid.
#'
#' @seealso This function is an updated variant of
#'   `CHBUtils::volcano_density_plot()`.
#'
#' @return Volcano plot arranged as grid (`grid = TRUE`), or [show()]
#'   individual [ggplot] (`grid = FALSE`).
#'
#' @examples
#' data(bcb)
#' dds <- DESeqDataSetFromTximport(
#'     txi = txi(bcb),
#'     colData = colData(bcb),
#'     design = formula(~group)) %>%
#'     DESeq
#' res <- results(dds)
#'
#' plot_volcano(res)
#' plot_volcano(res, genes = "Sulf1")
#' plot_volcano(res, ntop = 5L)
#' plot_volcano(res, padj = FALSE, alpha = 0.01, lfc = 4L)




#' @rdname plot_volcano
#' @param res_df [DESeqResults] coerced to a [data.frame]. By inputting a
#'   [data.frame], we allow for increased flexibility, since a user may want to
#'   subset the results prior to plotting.
.plot_volcano <- function(
    res_df,
    alpha = 0.1,
    padj = TRUE,
    lfc = 1L,
    genes = NULL,
    ntop = 0L,
    direction = "both",
    shade_color = "green",
    shade_alpha = 0.25,
    point_color = "gray",
    point_alpha = 0.75,
    point_outline_color = "darkgray",
    grid = TRUE) {
    if (!any(direction %in% c("both", "down", "up")) |
        length(direction) > 1L) {
        stop("Direction must be both, up, or down")
    }


    # Generate stats tibble ====
    stats <- res_df %>%
        rownames_to_column("ensgene") %>%
        as("tibble") %>%
        snake %>%
        tidy_select(!!!syms(c("ensgene",
                              "log2_fold_change",
                              "pvalue",
                              "padj"))) %>%
        # Filter zero counts for quicker plotting
        filter(!is.na(.data[["log2_fold_change"]]))
    g2s <- detectOrganism(stats[["ensgene"]][[1L]]) %>%
        gene2symbol
    stats <- left_join(stats, g2s, by = "ensgene")

    # Negative log10 transform the P values
    # Add `1e-10` here to prevent `Inf` values resulting from `log10()`
    if (isTRUE(padj)) {
        stats <- stats %>%
            mutate(neg_log10_pvalue = -log10(.data[["padj"]] + 1e-10))
        p_title <- "padj value"
    } else {
        stats <- stats %>%
            mutate(neg_log10_pvalue = -log10(.data[["pvalue"]] + 1e-10))
        p_title <- "p value"
    }

    stats <- stats %>%
        # Calculate rank score
        mutate(rank_score = .data[["neg_log10_pvalue"]] *
                   abs(.data[["log2_fold_change"]])) %>%
        arrange(desc(!!sym("rank_score")))


    # Text labels ====
    if (!is.null(genes)) {
        volcano_text <- stats %>%
            filter(.data[["symbol"]] %in% !!genes)
    } else if (ntop > 0L) {
        volcano_text <- stats[1L:ntop, ]
    } else {
        volcano_text <- NULL
    }


    # Plot ranges ====
    # Get range of LFC and P values to set up plot borders
    range_lfc <-
        c(floor(min(na.omit(stats[["log2_fold_change"]]))),
          ceiling(max(na.omit(stats[["log2_fold_change"]]))))
    range_neg_log10_pvalue <-
        c(floor(min(na.omit(stats[["neg_log10_pvalue"]]))),
          ceiling(max(na.omit(stats[["neg_log10_pvalue"]]))))


    # LFC density histogram ====
    lfc_density <- stats[["log2_fold_change"]] %>%
        na.omit %>%
        density
    lfc_density_df <- data.frame(x = lfc_density[["x"]],
                                 y = lfc_density[["y"]])
    lfc_hist <- stats %>%
        ggplot(aes_(x = ~log2_fold_change)) +
        geom_density() +
        scale_x_continuous(limits = range_lfc) +
        labs(x = "log2 fold change",
             y = "") +
        # Don't label density y-axis
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    if (direction == "both" | direction == "up") {
        lfc_hist <- lfc_hist +
            geom_ribbon(
                data = filter(lfc_density_df, .data[["x"]] > !!lfc),
                aes_(x = ~x, ymax = ~y),
                ymin = 0L,
                fill = shade_color,
                alpha = shade_alpha)
    }
    if (direction == "both" | direction == "down") {
        lfc_hist <- lfc_hist +
            geom_ribbon(
                data = filter(lfc_density_df, .data[["x"]] < -UQ(lfc)),
                aes_(x = ~x, ymax = ~y),
                ymin = 0L,
                fill = shade_color,
                alpha = shade_alpha)
    }


    # P value density plot ====
    pvalue_density <- stats[["neg_log10_pvalue"]] %>%
        na.omit %>%
        density
    pvalue_density_df <-
        data.frame(x = pvalue_density[["x"]],
                   y = pvalue_density[["y"]])
    pvalue_hist <- stats %>%
        ggplot(aes_(x = ~neg_log10_pvalue)) +
        geom_density() +
        geom_ribbon(data = filter(pvalue_density_df,
                                  .data[["x"]] > -log10(UQ(alpha) + 1e-10)),
                    aes_(x = ~x, ymax = ~y),
                    ymin = 0L,
                    fill = shade_color,
                    alpha = shade_alpha) +
        labs(x = paste("-log10", p_title),
             y = "") +
        # Don't label density y-axis
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())


    # Volcano plot ====
    volcano <- stats %>%
        ggplot(aes_(x = ~log2_fold_change,
                    y = ~neg_log10_pvalue)) +
        labs(x = "log2 fold change",
             y = paste("-log10", p_title)) +
        geom_point(
            alpha = point_alpha,
            color = point_outline_color,
            fill = point_color,
            pch = 21L) +
        theme(legend.position = "none") +
        scale_x_continuous(limits = range_lfc)
    if (!is.null(volcano_text)) {
        volcano <- volcano +
            geom_text_repel(
                data = volcano_text,
                aes_(x = ~log2_fold_change,
                     y = ~neg_log10_pvalue,
                     label = ~symbol),
                arrow = arrow(length = unit(0.01, "npc")),
                box.padding = unit(0.5, "lines"),
                color = "black",
                fontface = "bold",
                force = 1L,
                point.padding = unit(0.75, "lines"),
                segment.color = "gray",
                segment.size = 0.5,
                show.legend = FALSE,
                size = 4L)
    }
    if (direction == "both" | direction == "up") {
        volcano_poly_up <- with(stats, data.frame(
            x = as.numeric(c(
                lfc,
                lfc,
                max(range_lfc),
                max(range_lfc))),
            y = as.numeric(c(
                -log10(alpha + 1e-10),
                max(range_neg_log10_pvalue),
                max(range_neg_log10_pvalue),
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
                max(range_neg_log10_pvalue),
                max(range_neg_log10_pvalue),
                -log10(alpha + 1e-10)))))
        volcano <- volcano +
            geom_polygon(
                data = volcano_poly_down,
                aes_(x = ~x, y = ~y),
                fill = shade_color,
                alpha = shade_alpha)
    }


    # Grid layout ====
    if (isTRUE(grid)) {
        ggdraw() +
            # Coordinates are relative to lower left corner
            draw_plot(
                lfc_hist,
                x = 0L, y = 0.7, width = 0.5, height = 0.3) +
            draw_plot(
                pvalue_hist,
                x = 0.5, y = 0.7, width = 0.5, height = 0.3) +
            draw_plot(
                volcano, x = 0L, y = 0L, width = 1L, height = 0.7)
    } else {
        show(lfc_hist)
        show(pvalue_hist)
        show(volcano)
    }
}



#' @rdname plot_volcano
#' @export
setMethod("plot_volcano", "DESeqResults", function(object, alpha = NULL, ...) {
    res_df <- as.data.frame(object)
    if (is.null(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    .plot_volcano(res_df, alpha = alpha, ...)
})
