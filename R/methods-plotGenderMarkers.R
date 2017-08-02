#' Plot Sexually Dimorphic Gender Markers
#'
#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#'
#' @return [ggplot].
#' @export
#'
#' @examples
#' data(bcb)
#' plotGenderMarkers(bcb)
NULL



# Methods ====
#' @rdname plotGenderMarkers
#' @export
setMethod("plotGenderMarkers", "bcbioRNADataSet", function(object) {
    # Organism-specific dimorphic markers ====
    organism <- metadata(object)[["organism"]]

    # Load the relevant internal gender markers data
    if (organism == "mmusculus") {
        genderMarkers <- get("genderMarkersMmusculus",
                              envir = as.environment("package:bcbioRnaseq"))
    } else if (organism == "hsapiens") {
        stop("Draft support coming soon")
    } else {
        stop("Unsupported organism")
    }

    # Prepare the source gender markers data frame
    genderMarkers <- genderMarkers %>%
        .[.[["include"]] == TRUE, ]

    # Ensembl identifiers
    ensgene <- genderMarkers %>%
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

    tpm %>%
        as("tibble") %>%
        # Can also declare `measure.vars` here
        # If you don't set `id`, function will output a message
        melt(id = 1L) %>%
        setNames(c("ensgene", "sampleName", "counts")) %>%
        left_join(genderMarkers, by = "ensgene") %>%
        ggplot(aes_(x = ~symbol,
                    y = ~counts,
                    color = ~sampleName,
                    shape = ~chromosome)) +
        geom_jitter(size = 4L) +
        expand_limits(y = 0L) +
        labs(title = "gender markers",
             x = "gene",
             y = "transcripts per million (tpm)")
})
