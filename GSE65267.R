# GSE65267
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65267

# v0.2.4
# BiocInstaller::biocLite("hbc/bcbioRNASeq")
library(bcbioRNASeq)

# bcbio RNA-seq run is saved on HMS O2 cluster
gse65267 <- loadRNASeq(
    uploadDir = file.path(
        "/n"
        "data1",
        "cores",
        "bcbio",
        "bcbioRNASeq",
        "F1000v2",
        "GSE65267-merged",
        "final"
    ),
    organism = "Mus musculus",
    ensemblVersion = 90L
)

# Subset GSE65267 to day 0, 1, 3, 7 samples
bcb <- selectSamples(gse65267, day = c(0L, 1L, 3L, 7L))

saveData(gse65267, bcb, dir = "data")
