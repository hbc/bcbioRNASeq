library(bcbioRNASeq)
bcb <- loadRNASeq(
    uploadDir = "bcbio_rnaseq_run/final",
    interestingGroups = c("genotype", "treatment"),
    organism = "Homo sapiens"
)
saveData(bcb, dir = file.path("data", Sys.Date()))
