# GSE65267
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65267

library(bcbioRNASeq)

# Import bcbio RNA-seq run data from HMS O2 cluster
bcb <- loadRNASeq(
    uploadDir = file.path(
        "/n"
        "data1",
        "cores",
        "bcbio",
        "bcbioRNASeq",
        "F1000v2",
        "GSE65267-merged",
        "final"),
    interestingGroups = "treatment",
    ensemblVersion = 90L
)
saveData(bcb)
