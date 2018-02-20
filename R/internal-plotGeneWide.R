#' @importFrom dplyr arrange group_by left_join
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom magrittr set_colnames
#' @importFrom rlang !! !!! sym syms
#' @importFrom reshape2 melt
#' @importFrom tibble as_tibble rownames_to_column
.plotGeneWide <- function(
    object,
    genes,
    gene2symbol = NULL,
    colData,
    interestingGroups = "sampleName",
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = NULL) {
    assert_is_matrix(object)
    assert_has_dimnames(object)
    assert_is_character(genes)
    assert_is_subset(genes, rownames(object))
    object <- object[genes, , drop = FALSE]
    assertFormalGene2symbol(object, genes, gene2symbol)
    assertFormalAnnotationCol(object, colData)
    assertFormalInterestingGroups(colData, interestingGroups)
    assert_is_a_string(countsAxisLabel)
    assert_is_a_bool(medianLine)
    assertIsColorScaleDiscreteOrNULL(color)
    assertIsAStringOrNULL(title)

    colData <- uniteInterestingGroups(colData, interestingGroups)

    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        match <- match(x = genes, table = gene2symbol[["ensgene"]])
        symbols <- gene2symbol[match, "symbol", drop = TRUE]
        rownames(object) <- symbols
    }

    # Melt counts into long format
    data <- object %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1L) %>%
        as_tibble() %>%
        set_colnames(c("gene", "sampleID", "counts")) %>%
        arrange(!!!syms(c("gene", "sampleID"))) %>%
        group_by(!!sym("gene")) %>%
        left_join(colData, by = "sampleID")

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "gene",
            y = "counts",
            color = "interestingGroups")
    ) +
        geom_point(size = 4L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
        labs(
            title = title,
            x = "gene",
            y = countsAxisLabel,
            color = paste(interestingGroups, collapse = ":\n")
        ) +
        expand_limits(y = 0L)

    if (isTRUE(medianLine)) {
        p <- p + geneMedianLine
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(color = FALSE)
    }

    p
}
