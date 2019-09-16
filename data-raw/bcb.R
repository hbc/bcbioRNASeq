## bcbioRNASeq example objects.
## Updated 2019-09-16.

library(usethis)
library(pryr)

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

object_size(bcb)
## 889 kB

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
