# bcbioRNASeq example objects.
# Updated 2019-07-18.

library(pryr)
library(usethis)

# Restrict to 2 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(2e6, class = "object_size")

# GSE65267 =====================================================================
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

# F1000 paper ==================================================================
# Subset days 0, 1, 3, 7.
f1000 <- selectSamples(gse65267, day = c(0L, 1L, 3L, 7L))
saveData(f1000, dir = "~")

# bcb ====================================================================
# Minimal working example: days 0, 7.
bcb <- selectSamples(gse65267, day = c(0L, 7L))

# Minimize metadata slots that are too large for a working example.
metadata(bcb)[["bcbioCommandsLog"]] <- character()
metadata(bcb)[["bcbioLog"]] <- character()
metadata(bcb)[["dataVersions"]] <- DataFrame()
metadata(bcb)[["programVersions"]] <- DataFrame()
metadata(bcb)[["tx2gene"]] <- head(metadata(bcb)[["tx2gene"]])
metadata(bcb)[["yaml"]] <- list()

# Always include sexually dimorphic marker genes.
dimorphic_genes <- c(
    "ENSMUSG00000056673",
    "ENSMUSG00000068457",
    "ENSMUSG00000069045",
    "ENSMUSG00000086503"
)
stopifnot(all(dimorphic_genes %in% rownames(bcb)))
# Include other genes ranked by abundance.
top_genes <- tpm(bcb) %>%
    rowSums() %>%
    sort(decreasing = TRUE) %>%
    names()
genes <- c(dimorphic_genes, top_genes) %>%
    unique() %>%
    head(n = 500L) %>%
    sort()
bcb <- bcb[genes, ]

# Update the interesting groups and design formula.
interestingGroups(bcb) <- "treatment"

# DESeq2 doesn't like spaces in design formula factors.
bcb[["treatment"]] <- snake(bcb[["treatment"]])

# Report the size of each slot in bytes.
vapply(
    X = coerceS4ToList(bcb),
    FUN = object_size,
    FUN.VALUE = numeric(1L)
)
object_size(bcb)
stopifnot(object_size(bcb) < limit)

# Check that object is valid.
stopifnot(is(bcb, "bcbioRNASeq"))
validObject(bcb)

usethis::use_data(bcb, overwrite = TRUE, compress = "xz")
