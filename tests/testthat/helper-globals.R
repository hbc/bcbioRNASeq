## nolint start

## Fix for pheatmap partial match warning.
## https://github.com/raivokolde/pheatmap/issues/46
options(
    "warnPartialMatchAttr" = FALSE,
    "warnPartialMatchDollar" = FALSE
)

data(bcb, envir = environment())
## > invisible(validObject(bcb))

bcb_fast <- readRDS(file.path("cache", "bcb_fast.rds")) # nolint
## > invisible(validObject(bcb_fast))

Gene2Symbol <- AcidGenerics::Gene2Symbol
`sampleData<-` <- AcidGenerics::`sampleData<-`
assay <- SummarizedExperiment::assay
hasInternet <- goalie::hasInternet
initDir <- AcidBase::initDir
nonzeroRowsAndCols <- AcidGenerics::nonzeroRowsAndCols
render <- rmarkdown::render
seqnames <- GenomeInfoDb::seqnames
tempdir2 <- AcidBase::tempdir2

object <- bcb
g2s <- Gene2Symbol(object)
geneIds <- head(g2s[["geneId"]])
geneNames <- head(g2s[["geneName"]])

uploadDir <- system.file("extdata", "bcbio", package = "bcbioRNASeq")
organism <- "Mus musculus"
ensemblRelease <- 90L

## nolint end
