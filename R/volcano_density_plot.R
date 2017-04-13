# https://github.com/hbc/CHBUtils/blob/master/R/volcanoPlot.R



# volcano_density_plot: no visible global function definition for
# ‘density’
# volcano_density_plot: no visible binding for global variable ‘logFC’
# volcano_density_plot: no visible binding for global variable ‘V1’
# volcano_density_plot: no visible binding for global variable ‘V2’
# volcano_density_plot: no visible binding for global variable
# ‘adj.P.Val’
# volcano_density_plot: no visible binding for global variable ‘x’
# volcano_density_plot: no visible binding for global variable ‘y’
# volcano_density_plot: no visible binding for global variable ‘name’
# volcano_density_plot: no visible global function definition for
# ‘grid.arrange’
# Undefined global functions or variables:
#     V1 V2 adj.P.Val asis_output density grid.arrange import_metadata
# logFC name opts_knit x y



fmt <- function() {
    function(x) format(x, nsmall = 1, scientific = FALSE)
}



#' Create volcano plot from log2FC and adjusted pvalues data frame
#'
#' @param stats data.frame with two columns: logFC and Adjusted.Pvalue
#' @param side plot UP, DOWN or BOTH de-regulated points
#' @param title title for the figure
#' @param pval_cutoff cutoff for the adjusted pvalue. Default 0.05
#' @param lfc_cutoff cutoff for the log2FC. Default 1
#' @param shade_alpha transparency value. Default 0.25
#' @param point_color colours for points. Default gray
#' @param point_alpha transparency for points. Default 0.75
#' @param point_outline_color Default darkgray
#' @param line_color Defaul gray
#' @param plot_text Data frame similar to stats table, with a third column being
#'   the text you want to add next to the point in the volcano plot.
#' @export
volcano_density_plot <- function(
    stats,
    side = "both",
    title = "volcano plot with marginal distributions",
    pval_cutoff = 0.05,
    lfc_cutoff = 1,
    shade.colour = "green",
    shade_alpha = 0.25,
    point_color = "gray",
    point_alpha = 0.75,
    point_outline_color = "darkgray",
    line_color = "gray",
    plot_text = NULL) {
    if (!any(side %in% c("both", "down", "up")) | length(side) > 1)
        stop("side parameter should be: both, up or down")
    if (ncol(stats) != 2)
        stop("need a data frame with two columns: logFC and Adjusted.Pvalue")
    if (sum(rowSums(is.na(stats))) > 0)
        stats = stats[rowSums(is.na(stats)) == 0, ]
    if (any(stats[, 2] > 1) | any(stats[, 2] < 0) )
        stop("pvalues needs to be > 0 and < 1")
    names(stats) <- c("logFC", "adj.P.Val")
    stats[, 2] <- stats[, 2] + 1e-10

    # Get range of log fold change and p-value values to set up plot borders
    range.lfc <- c(floor(min(stats$logFC)) - 1,
                   ceiling(max(stats$logFC)) + 1)
    range.pval <- c(floor(min(-log10(stats$adj.P.Val))),
                    ceiling(max(-log10(stats$adj.P.Val))) + 1)

    # Make top plot: density plot with fold changes
    lfcd <- as.data.frame(cbind(density(stats$logFC)$x,
                                density(stats$logFC)$y))
    hist_top <- ggplot(data = stats, aes(x = logFC)) +
        geom_density(color = line_color) +
        theme(axis.title.x = element_blank()) +
        theme(plot.margin = unit(c(3, -5.5, 4, 3), "mm")) +
        scale_x_continuous(limits = range.lfc,
                           breaks = range.lfc[1]:range.lfc[2],
                           expand = c(.05,.05)) +
        scale_y_continuous(labels = fmt()) +
        labs(list(title = "fold changes density plot"))

    if (side == "both" | side == "up") {
        hist_top <- hist_top +
            geom_ribbon(data = subset(lfcd, V1 > lfc_cutoff),
                        aes(x = V1,ymax = V2),
                        ymin = 0,
                        fill = shade.colour,
                        alpha = shade_alpha)
    }
    if (side == "both" | side == "down") {
        hist_top <- hist_top +
            geom_ribbon(data = subset(lfcd, V1 < -lfc_cutoff),
                        aes(x = V1, ymax = V2),
                        ymin = 0,
                        fill = shade.colour,
                        alpha = shade_alpha)
    }

    # Make blank plot
    empty <- ggplot() +
        geom_point(aes(1, 1), colour = "white") +
        theme(panel.grid = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()
        )

    # Make scatter volcano plot
    scat.poly.up <- with(stats, data.frame(
        x = as.numeric(c(lfc_cutoff,
                         lfc_cutoff,
                         max(range.lfc),
                         max(range.lfc))),
        y = as.numeric(c(-log10(pval_cutoff),
                         max(range.pval),
                         max(range.pval),
                         -log10(pval_cutoff)))))
    scat.poly.down <- with(stats, data.frame(
        x = as.numeric(c(-lfc_cutoff,
                         -lfc_cutoff,
                         min(range.lfc),
                         min(range.lfc))),
        y = as.numeric(c(-log10(pval_cutoff),
                         max(range.pval),
                         max(range.pval),
                         -log10(pval_cutoff)))))
    scatter <- ggplot(data = stats,
                      aes(x = logFC,
                          y = -log10(adj.P.Val))) +
        geom_point(alpha = point_alpha,
                   pch = 21,
                   fill = point_color,
                   color = point_outline_color) +
        xlab("log2 fold change") +
        ylab("-log10(adjusted p-value)") +
        theme(legend.position = "none") +
        theme(plot.margin = unit(c(3,-5.5,4,3), "mm") ) +
        scale_x_continuous(limits = range.lfc,
                           breaks = range.lfc[1]:range.lfc[2],
                           expand = c(.05,.05)) +
        scale_y_continuous(labels = fmt(),
                           limits = range.pval) +
        labs(list(title = "Volcano plot"))

    if (side == "both" | side == "up") {
        scatter <- scatter +
            geom_polygon(data = scat.poly.up,
                         aes(x = x,
                             y = y),
                         fill = shade.colour,
                         alpha = shade_alpha)
    }

    if (side == "both" | side == "down") {
        scatter <- scatter +
            geom_polygon(data = scat.poly.down,
                         aes(x = x,
                             y = y),
                         fill = shade.colour,
                         alpha = shade_alpha)
    }

    if (!is.null(plot_text)) {
        names(plot_text) <- c("logFC", "adj.P.Val", "name")
        plot_text[, 2] <- plot_text[, 2] + 1e-10
        scatter <- scatter +
            geom_text_repel(data = plot_text,
                            aes(x = logFC,
                                y = -log10(adj.P.Val),
                                label = name),
                            size = 3)
    }

    # make right plot - density plot of adjusted pvalues
    pvald <- as.data.frame(
        cbind(density(-log10(stats$adj.P.Val))$x,
              density(-log10(stats$adj.P.Val))$y))
    hist_right <- ggplot(data = stats,
                         aes(x = -log10(adj.P.Val))) +
        geom_density(color = line_color) +
        geom_ribbon(data = subset(pvald, V1 > -log10(pval_cutoff)),
                    aes(x = V1,
                        ymax = V2),
                    ymin = 0,
                    fill = shade.colour,
                    alpha = shade_alpha) +
        coord_flip() +
        scale_x_continuous(limits = range.pval) +
        theme(axis.title.y = element_blank()) +
        theme(plot.margin = unit(c(3, -5.5, 4,3), "mm")) +
        labs(list(title = "adj.pval density plot"))

    # plot all plots
    pp.logfc <- ggplotGrob(hist_top)
    pp.empty <- ggplotGrob(empty)
    pp.volc <- ggplotGrob(scatter)
    pp.pval <- ggplotGrob(hist_right)

    plot <- grid.arrange(
        top = textGrob(title),
        arrangeGrob(pp.logfc, pp.volc,
                    heights = c(1, 3),
                    col = 1),
        arrangeGrob(pp.empty, pp.pval,
                    heights = c(1, 3),
                    ncol = 1),
        ncol = 2,
        widths = c(3, 1))

    return(plot)
}
