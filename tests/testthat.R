set.seed(42L)

library(testthat)
library(bcbioRNASeq)
data(bcb, dds, res)

test_check("bcbioRNASeq")
