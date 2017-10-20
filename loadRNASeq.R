# GSE65267 example dataset
# Run this script on O2
library(bcbioRNASeq)
bcb <- loadRNASeq(
    uploadDir = file.path(
        "/n",
        "scratch2",
        "hsph_bioinformatic_core",
        "lp113",
        "workflow",
        "samples-merged",
        "final"),
    interestingGroups = "group",
    ensemblVersion = "current"
)
saveData(bcb)
