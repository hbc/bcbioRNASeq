#' @rdname tmm
.tmm <- function(object) {
    message("Generating TMM-normalized counts with edgeR")
    object %>%
        as.matrix %>%
        DGEList %>%
        calcNormFactors %>%
        cpm(normalized.lib.sizes = TRUE)
}
