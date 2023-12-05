## nolint start

## Fix for pheatmap partial match warning.
## https://github.com/raivokolde/pheatmap/issues/46
options(
    "warnPartialMatchAttr" = FALSE,
    "warnPartialMatchDollar" = FALSE
)

GeneToSymbol <- AcidGenerics::GeneToSymbol
`design<-` <- BiocGenerics::`design<-`
`sampleData<-` <- AcidGenerics::`sampleData<-`
assay <- SummarizedExperiment::assay
data <- utils::data
hasInternet <- goalie::hasInternet
initDir <- AcidBase::initDir
isADir <- goalie::isADir
nonzeroRowsAndCols <- AcidGenerics::nonzeroRowsAndCols
pasteUrl <- AcidBase::pasteUrl
seqnames <- GenomeInfoDb::seqnames
untar <- utils::untar

data(bcb, envir = environment())
bcb_fast <- readRDS(file.path(cacheDir, "bcb_fast.rds")) # nolint

object <- bcb
g2s <- GeneToSymbol(object)
geneIds <- head(g2s[["geneId"]])
geneNames <- head(g2s[["geneName"]])

uploadDir <- system.file("extdata", "bcbio", package = "bcbioRNASeq")
organism <- "Mus musculus"
ensemblRelease <- 90L

## nolint end
