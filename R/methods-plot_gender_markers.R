#' Plot sexually dimorphic gender markers
#'
#' @rdname plot_gender_markers
#' @author Michael Steinbaugh
#'
#' @param object Object.
#'
#' @export
setMethod("plot_gender_markers", "bcbioRNADataSet", function(object) {
    # Organism-specific dimorphic markers ====
    organism <- metadata(object)[["organism"]]

    # Load the relevant internal gender markers data
    if (organism == "mmusculus") {
        gender_markers <- get("gender_markers_mmusculus",
                              envir = as.environment("package:bcbioRnaseq"))
    } else {
        # FIXME Need to add support for human markers
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


    # ggplot ====
    tpm(object) %>%
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
})
