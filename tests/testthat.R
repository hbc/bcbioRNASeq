library(testthat)
library(bcbioRnaseq)

set.seed(42)

data(bcb)
dds <- DESeqDataSetFromTximport(
    txi = txi(bcb),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq
res <- results(dds)

test_check("bcbioRnaseq")
