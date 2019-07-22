## bcbioRNASeq example objects.
## Updated 2019-07-19.

library(pryr)
library(usethis)

## Restrict to 2 MB.
## Use `pryr::object_size()` instead of `utils::object.size()`.
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

## FIXME The levels of the rowRanegs need to be resized in makeSummarizedExperiment call.

object_size(bcb)
## 11.1 MB

## The rowRanges are quite large. Need to resize the levels automatically.

lapply(coerceS4ToList(bcb), object_size)

## Check the size of the metadata.
lapply(metadata(bcb), object_size)

## Minimize metadata slots that are too large for a working example.
metadata(bcb)[["bcbioCommandsLog"]] <- character()
metadata(bcb)[["bcbioLog"]] <- character()
metadata(bcb)[["dataVersions"]] <- DataFrame()
metadata(bcb)[["programVersions"]] <- DataFrame()

metadata(bcb)[["yaml"]] <- list()

## Update the interesting groups and design formula.
interestingGroups(bcb) <- "treatment"

## DESeq2 doesn't like spaces in design formula factors.
bcb[["treatment"]] <- snake(bcb[["treatment"]])

## Report the size of each slot in bytes.
lapply(coerceS4ToList(bcb), object_size)
object_size(bcb)
stopifnot(object_size(bcb) < limit)

## Check that object is valid.
stopifnot(is(bcb, "bcbioRNASeq"))
validObject(bcb)

usethis::use_data(bcb, overwrite = TRUE, compress = "xz")
