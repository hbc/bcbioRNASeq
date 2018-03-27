library(pryr)
library(devtools)
library(tidyverse)
load_all()

# GSE65267
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65267
# HMS O2: /n/data1/cores/bcbio/bcbioRNASeq/F1000v2/GSE65267-merged/final
gse65267 <- loadRNASeq(
    uploadDir = file.path(
        "~",
        "O2",
        "bcbio",
        "bcbioRNASeq",
        "F1000v2",
        "GSE65267-merged",
        "final"
    ),
    level = "genes",
    caller = "salmon",
    interestingGroups = c("treatment", "day"),
    organism = "Mus musculus",
    genomeBuild = "GRCm38",
    ensemblRelease = 90L
)
# Too large to save in the package
saveData(gse65267, dir = "~")

# F1000 paper: days 0, 1, 3, 7
# Working example: days 0, 7
bcb_small <- selectSamples(gse65267, day = c(0L, 7L), transform = FALSE)
stopifnot(identical(
    names(assays(bcb_small)),
    c("raw", "tpm", "length")
))

# Minimize metadata slots that take up disk space
metadata(bcb_small)[["bcbioCommandsLog"]] <- character()
metadata(bcb_small)[["bcbioLog"]] <- character()
metadata(bcb_small)[["dataVersions"]] <- tibble()
metadata(bcb_small)[["programVersions"]] <- tibble()
metadata(bcb_small)[["tx2gene"]] <- head(metadata(bcb_small)[["tx2gene"]])
metadata(bcb_small)[["yaml"]] <- list()

# Sort the genes by abundance
abundance <- tpm(bcb_small) %>%
    rowSums() %>%
    sort(decreasing = TRUE) %>%
    names()
# Ensure all dimorphic gender markers are included
dimorphic <- genderMarkers[["musMusculus"]] %>%
    filter(include == TRUE) %>%
    pull(geneID) %>%
    sort() %>%
    intersect(rownames(bcb_small))
# Subset to include only the top genes of interest
genes <- c(dimorphic, abundance) %>%
    unique() %>%
    head(500L)
bcb_small <- bcb_small[genes, , transform = TRUE]
stopifnot(identical(
    names(assays(bcb_small)),
    c("raw", "tpm", "length", "rlog", "vst")
))

# Update the interesting groups and design formula
interestingGroups(bcb_small) <- "treatment"

# DESeq2 doesn't like spaces in design factors, so fix that in colData
bcb_small[["treatment"]] <- snake(bcb_small[["treatment"]])

# Check that object is valid
validObject(bcb_small)

# Check the object size
pryr::object_size(bcb_small)
format(object.size(bcb_small), units = "auto")

# Check sizes in slots
lapply(assays(bcb_small), pryr::object_size)
lapply(metadata(bcb_small), pryr::object_size)

lapply(assays(bcb_small), object.size)
lapply(metadata(bcb_small), object.size)

use_data(bcb_small, overwrite = TRUE, compress = "xz")
