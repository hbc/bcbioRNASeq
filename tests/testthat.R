set.seed(42L)

library(testthat)
library(bcbioRnaseq)
data(bcb, dds, res)

test_check("bcbioRnaseq")
