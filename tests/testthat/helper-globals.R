## Fix for pheatmap partial match warning.
## https://github.com/raivokolde/pheatmap/issues/46
options(
    warnPartialMatchAttr = FALSE,
    warnPartialMatchDollar = FALSE
)

data(bcb, envir = environment())
invisible(validObject(bcb))

object <- bcb
g2s <- Gene2Symbol(object)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
organism <- "Mus musculus"
ensemblRelease <- 90L

## nolint start
assay <- SummarizedExperiment::assay
hasInternet <- goalie::hasInternet
initDir <- basejump::initDir
skip_on_docker <- goalie::skip_on_docker
## nolint end
