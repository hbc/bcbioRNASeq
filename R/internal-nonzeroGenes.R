.nonzeroGenes <- function(object) {
    raw <- counts(object, normalized = FALSE)
    nonzero <- rowSums(raw) > 0L
    genes <- rownames(raw[nonzero, , drop = FALSE])
    inform(paste(length(genes), "non-zero genes detected"))
    sort(genes)
}
