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

# Use day 0 and 7 in minimal example
# F1000 paper: days 0, 1, 3, 7
bcb_small <- selectSamples(gse65267, day = c(0L, 7L))

# Minimize metadata slots that take up disk space
metadata(bcb_small)$bcbioCommandsLog <- character()
metadata(bcb_small)$bcbioLog <- character()
metadata(bcb_small)$dataVersions <- tibble()
metadata(bcb_small)$programVersions <- tibble()
metadata(bcb_small)$tx2gene <- head(metadata(bcb_small)$tx2gene)
metadata(bcb_small)$yaml <- list()

# Sort the genes by abundance
abundance <- tpm(bcb_small) %>%
    rowSums() %>%
    sort(decreasing = TRUE) %>%
    names()
# Ensure all dimorphic gender markers are included
dimorphic <- genderMarkers$musMusculus %>%
    filter(include == TRUE) %>%
    pull(geneID) %>%
    sort() %>%
    intersect(rownames(bcb_small))
# Subset to include 100 genes
genes <- c(dimorphic, abundance) %>%
    unique() %>%
    head(100L)
bcb_small <- bcb_small[genes, ]

# DESeq2 doesn't like spaces in design factors, so fix that first
bcb_small$treatment <- snake(bcb_small$treatment)
validObject(bcb_small)

# Check the object size
object.size(bcb_small) %>% format(units = "auto")
pryr::object_size(bcb_small)

# Update the interesting groups and design formula
interestingGroups(bcb_small) <- "treatment"
design(bcb_small) <- formula(~treatment)
validObject(bcb_small)

# Check sizes
lapply(assays(bcb_small), object.size)
lapply(metadata(bcb_small), object.size)

dds_small <- assays(bcb_small)[["dds"]]
stopifnot(design(dds_small) == formula(~treatment))
stopifnot(identical(
    resultsNames(dds_small),
    c("Intercept", "treatment_folic_acid_vs_control")
))

# Include an example where metadata has a space, sanitized with `make.names()`
res_small <- results(
    dds_small,
    contrast = c(
        factor = "treatment",
        numerator = "folic_acid",
        denominator = "control"
    )
)

use_data(bcb_small, res_small, overwrite = TRUE, compress = "xz")
