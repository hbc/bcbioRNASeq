library(bcbioRNASeq)
bcb <- loadRNASeq(
    file.path("bcbio_rnaseq_run", "final"),
    interestingGroups = c("genotype", "treatment"))
# Back up all data inside bcbioRNASeq object
flatFiles <- flatFiles(bcb)
saveData(bcb, flatFiles, dir = "data")
