#' Plot Sexually Dimorphic Gender Markers
#'
#' @rdname plot_gender_markers
#' @author Michael Steinbaugh
#'
#' @return [ggplot].
#' @export
#'
#' @examples
#' data(bcb)
#' plot_gender_markers(bcb)
setMethod("plot_gender_markers", "bcbioRNADataSet", function(object) {
    # Organism-specific dimorphic markers ====
    organism <- metadata(object)[["organism"]]

    # Load the relevant internal gender markers data
    if (organism == "mmusculus") {
        gender_markers <- get("gender_markers_mmusculus",
                              envir = as.environment("package:bcbioRnaseq"))
    } else if (organism == "hsapiens") {
        stop("Draft support coming soon")
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

    # Subset TPM
    tpm <- tpm(object)
    if (!all(ensgene %in% rownames(tpm))) {
        warning("Missing dimorphic genes in counts matrix")
        return(NULL)
    }
    tpm <- tpm[ensgene, ]

    # Return ggplot
    tpm %>%
        as("tibble") %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1L) %>%
        set_names(c("ensgene", "sample_name", "counts")) %>%
        left_join(gender_markers, by = "ensgene") %>%
        ggplot(aes_(x = ~symbol,
                    y = ~counts,
                    color = ~sample_name,
                    shape = ~chromosome)) +
        geom_jitter(size = 4L) +
        expand_limits(y = 0L) +
        labs(title = "gender markers",
             x = "gene",
             y = "transcripts per million (tpm)")
})
