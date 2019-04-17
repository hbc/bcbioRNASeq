data(bcb, envir = environment())

object <- bcb
g2s <- Gene2Symbol(object)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])

# nolint start
assay <- SummarizedExperiment::assay
# nolint end
