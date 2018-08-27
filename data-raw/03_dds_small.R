#' DESeq2 example objects
#' Last updated 2018-08-27

library(DESeq2)

# DESeqDataSet =================================================================
dds_small <- as(bcb_small, "DESeqDataSet")
design(dds_small) <- ~ treatment
dds_small <- DESeq(dds_small)
validObject(dds_small)

stopifnot(design(dds_small) == ~ treatment)
stopifnot(identical(
    resultsNames(dds_small),
    c("Intercept", "treatment_folic_acid_vs_control")
))

# DESeqResults =================================================================
resultsNames(dds_small)
res_small <- results(
    dds_small,
    contrast = c(
        factor = "treatment",
        numerator = "folic_acid",
        denominator = "control"
    )
)
# Shrink log2 fold changes
res_small <- lfcShrink(
    dds = dds_small,
    coef = 2L,
    res = res_small
)
validObject(res_small)

# Save =========================================================================
devtools::use_data(
    dds_small, res_small,
    overwrite = TRUE, compress = "xz"
)
