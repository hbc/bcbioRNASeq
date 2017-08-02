dds <- DESeqDataSetFromTximport(
    txi = txi(bcb),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq
saveData(dds, compress = "xz")
