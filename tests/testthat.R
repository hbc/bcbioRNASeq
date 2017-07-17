set.seed(42L)

library(testthat)
library(bcbioRnaseq)

data(bcb)
dds <- DESeqDataSetFromTximport(
    txi = txi(bcb),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq
res <- results(dds)

test_check("bcbioRnaseq")
