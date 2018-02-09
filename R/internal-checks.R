.checkGenes <- function(genes, gene2symbol) {
    if (!is.character(genes)) {
        abort("`genes` must be a character vector")
    }
    if (!(is.data.frame(gene2symbol) || is.null(gene2symbol))) {
        abort("`gene2symbol` must be a data.frame or NULL")
    }
    if (is.data.frame(gene2symbol)) {
        checkGene2symbol(gene2symbol)
        if (!all(genes %in% gene2symbol[["ensgene"]])) {
            abort(paste(
                "Genes missing in gene2symbol:",
                setdiff(genes, gene2symbol[["ensgene"]])
            ))
        }
    }
}
