library(bcbioRNASeq)
bcb <- loadRNASeqRun(
    file.path("bcbio_rnaseq_run", "final"),
    interestingGroups = c("genotype", "treatment"))
saveData(bcb)
