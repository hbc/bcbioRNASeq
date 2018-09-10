#' Example bcbioRNASeq object
#' Last updated 2018-09-09

library("tidyverse")

# FIXME Improve the object size limit method here (see basejump example)

# GSE65267 =====================================================================
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65267
# HMS O2: /n/data1/cores/bcbio/bcbioRNASeq/F1000v2
# Using sshfs connection to O2 here
gse65267 <- bcbioRNASeq(
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
    ensemblRelease = 90L,
    vst = TRUE,
    rlog = FALSE
)
saveData(gse65267, dir = "data-raw")

# F1000 paper ==================================================================
# Subset days 0, 1, 3, 7
f1000 <- selectSamples(gse65267, day = c(0L, 1L, 3L, 7L))
saveData(f1000, dir = "data-raw")

# bcb_small ====================================================================
# Minimal working example: days 0, 7
bcb <- selectSamples(gse65267, day = c(0L, 7L))

# Minimize metadata slots that are too large for a working example
metadata(bcb)[["bcbioCommandsLog"]] <- character()
metadata(bcb)[["bcbioLog"]] <- character()
metadata(bcb)[["dataVersions"]] <- tibble()
metadata(bcb)[["programVersions"]] <- tibble()
metadata(bcb)[["tx2gene"]] <- head(metadata(bcb)[["tx2gene"]])
metadata(bcb)[["yaml"]] <- list()

# Sort the genes by abundance
abundance <- tpm(bcb) %>%
    rowSums() %>%
    sort(decreasing = TRUE) %>%
    names()
# Ensure all dimorphic gender markers are included
dimorphic <- bcbioRNASeq::gender_markers[["musMusculus"]] %>%
    pull(geneID) %>%
    sort() %>%
    intersect(rownames(bcb))
# Subset to include only the top genes of interest
genes <- c(dimorphic, abundance) %>%
    unique() %>%
    head(500L) %>%
    sort()
bcb <- bcb[genes, ]

# Update the interesting groups and design formula
interestingGroups(bcb) <- "treatment"

# DESeq2 doesn't like spaces in design factors, so fix that in colData
bcb[["treatment"]] <- snake(bcb[["treatment"]])

# Check that object is valid
stopifnot(is(bcb, "bcbioRNASeq"))
validObject(bcb)

# Check the object size
format(object.size(bcb), units = "auto")
pryr::object_size(bcb)
lapply(assays(bcb), pryr::object_size)
lapply(metadata(bcb), pryr::object_size)

bcb_small <- bcb
devtools::use_data(bcb_small, overwrite = TRUE, compress = "xz")
