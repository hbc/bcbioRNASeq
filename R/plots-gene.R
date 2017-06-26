# TODO Look at incorporating [plotCounts()] into this function.
# TODO Function needs rework for annotables
#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param bcb [bcbioRNADataSet].
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param format Ensembl identifier format. Defaults to the gene symbol (a.k.a.
#'   `external_gene_name`).
#'
#' @export
plot_gene <- function(
    bcb,
    gene,
    format = "symbol") {
    if (!format %in% c("ensgene", "symbol")) {
        stop("Unsupported gene identifier format")
    }
    counts <- tpm(bcb)
    metadata <- colData(bcb)
    color <- metadata(bcb)[["interesting_groups"]][[1L]]

    # Match unique gene identifier with name (gene symbol) using the internally
    # stored Ensembl annotations saved in the run object
    match <- metadata(bcb)[["annotable"]] %>%
        as.data.frame %>%
        filter(.data[[format]] %in% gene) %>%
        arrange(!!sym("symbol"))

    # Seq along Ensembl data frame here instead of the gene input vector,
    # which will then output only genes that match Ensembl
    lapply(1L:nrow(match), function(a) {
        gene_name <- match[["symbol"]][[a]]
        gene_id <- match[["ensgene"]][[a]]
        plot <- data.frame(
            x = colnames(counts),
            y = counts[gene_id, ],
            color = metadata[[color]]) %>%
            ggplot(
                aes_(x = ~x,
                     y = ~y,
                     color = ~color,
                     shape = ~color)
            ) +
            ggtitle(paste(gene_name,
                          paste0("(", gene_id, ")"))) +
            geom_point(size = 4L) +
            theme(
                axis.text.x = element_text(angle = 90L)
            ) +
            labs(x = "sample",
                 y = "transcripts per million (tpm)",
                 fill = "") +
            expand_limits(y = 0L)
        show(plot)
    }) %>% invisible
}



#' Plot sexually dimorphic gender markers
#'
#' @param bcb [bcbioRNADataSet].
#'
#' @export
plot_gender_markers <- function(bcb) {
    # Organism-specific dimorphic markers ----
    organism <- metadata(bcb)[["organism"]]

    # Load the relevant internal gender markers data
    if (organism == "mmusculus") {
        gender_markers <- get("gender_markers_musculus",
                              envir = as.environment("package:bcbioRnaseq"))
    } else {
        stop("Unsupported organism")
    }

    # Prepare the source gender markers data frame
    gender_markers <- gender_markers %>%
        filter(.data[["include"]] == TRUE)

    # Ensembl identifiers
    ensgene <- gender_markers %>%
        pull("ensgene") %>%
        sort %>%
        unique


    # ggplot ----
    tpm(bcb) %>%
        .[ensgene, ] %>%
        as.data.frame %>%
        rownames_to_column %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1L) %>%
        set_names(c("ensgene",
                    "description",
                    "counts")) %>%
        left_join(gender_markers, by = "ensgene") %>%
        ggplot(
            aes_(x = ~symbol,
                 y = ~counts,
                 color = ~description,
                 shape = ~chromosome)) +
        ggtitle("gender markers") +
        geom_jitter(size = 4L) +
        expand_limits(y = 0L) +
        labs(x = "gene",
             y = "transcripts per million (tpm)")
}
