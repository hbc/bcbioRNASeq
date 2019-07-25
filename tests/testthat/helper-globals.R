## nolint start
assay <- SummarizedExperiment::assay
skip_on_docker <- goalie::skip_on_docker
## nolint end

data(bcb, envir = environment())
invisible(validObject(bcb))

object <- bcb
g2s <- Gene2Symbol(object)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
organism <- "Mus musculus"
ensemblRelease <- 90L
