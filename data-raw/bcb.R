## bcbioRNASeq example object.
## Updated 2022-03-07.
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(lobstr)
    library(basejump)
})
load_all()
## Restrict to 2 MB.
limit <- structure(2e6, class = "object_size")
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
bcb <- bcbioRNASeq(
    uploadDir = uploadDir,
    level = "genes",
    caller = "salmon",
    interestingGroups = c("treatment", "day"),
    organism = "Mus musculus",
    ensemblRelease = 90L
)
## DESeq2 doesn't like spaces in design formula factors.
bcb[["treatment"]] <- snake(bcb[["treatment"]])
## Report the size of each slot in bytes.
lapply(coerceToList(bcb), obj_size)
stopifnot(
    isTRUE(obj_size(bcb) < limit),
    is(bcb, "bcbioRNASeq"),
    validObject(bcb)
)
usethis::use_data(bcb, overwrite = TRUE, compress = "xz")
