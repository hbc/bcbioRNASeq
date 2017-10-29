set.seed(42L)

library(testthat)
library(bcbioRNASeq)

# Use F1000 workflow example data (GSE65267) for unit tests
# https://github.com/hbc/bcbioRNASeq/tree/f1000v1
loadRemoteData(
    file.path(
        "https://github.com", "hbc", "bcbioRNASeq",
        "raw", "f1000v1", "data", "bcb.rda"),
    quiet = TRUE)
gse65267 <- bcb
rm(bcb)

test_check("bcbioRNASeq")
