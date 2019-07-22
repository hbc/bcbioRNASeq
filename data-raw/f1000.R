# F1000 workflow paper example.
# Updated 2019-07-19.

library(pryr)
library(usethis)

# Restrict to 2 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65267
# HMS O2: /n/data1/cores/bcbio/bcbioRNASeq/F1000v2
# Using sshfs connection to O2 here.
gse65267 <- bcbioRNASeq(
    uploadDir = file.path(
        "",
        "mnt",
        "O2",
        "hbc",
        "bcbioRNASeq",
        "F1000v2",
        "GSE65267-merged",
        "final"
    ),
    level = "genes",
    caller = "salmon",
    interestingGroups = c("treatment", "day"),
    organism = "Mus musculus",
    ensemblRelease = 90L
)
saveData(gse65267, dir = "~")

# Subset days 0, 1, 3, 7.
f1000 <- selectSamples(gse65267, day = c(0L, 1L, 3L, 7L))
saveData(f1000, dir = "~")
