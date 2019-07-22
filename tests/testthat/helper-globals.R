## nolint start
assay <- SummarizedExperiment::assay
skip_on_docker <- goalie::skip_on_docker
## nolint end

data(bcb, envir = environment())

object <- bcb
g2s <- Gene2Symbol(object)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])
