#' Example DESeq2 differential expression analysis
#' Last updated 2018-10-08

library(DESeq2)

# Restrict to 2 MB.
limit <- structure(2e6, class = "object_size")

# DESeqDataSet
# Coerce from bcbioRNASeq object.
dds <- as(bcb_small, "DESeqDataSet")
design(dds) <- ~ treatment
dds <- DESeq(dds)
validObject(dds)
stopifnot(design(dds) == ~ treatment)
stopifnot(identical(
    resultsNames(dds),
    c("Intercept", "treatment_folic_acid_vs_control")
))

# DESeqTransform
vst <- varianceStabilizingTransformation(dds)

# DESeqResults
resultsNames(dds)
res <- results(
    object = dds,
    contrast = c(
        factor = "treatment",
        numerator = "folic_acid",
        denominator = "control"
    )
)

# Shrink log2 fold changes
# `BiocManager::install("apeglm")`
res_shrunken <- lfcShrink(
    dds = dds,
    coef = 2L,
    res = res,
    type = "apeglm"
)
validObject(res_shrunken)

deseq_small <- DESeqAnalysis(
    data = dds,
    transform = vst,
    results = list(res),
    lfcShrink = list(res_shrunken)
)
validObject(deseq_small)
print(deseq_small)

# Report the size of each slot in bytes.
vapply(
    X = coerceS4ToList(deseq_small),
    FUN = object.size,
    FUN.VALUE = numeric(1L)
)
format(object.size(deseq_small), units = "auto")
stopifnot(object.size(bcb) < limit)

# Check that object is valid.
stopifnot(is(deseq_small, "DESeqAnalysis"))
stopifnot(validObject(deseq_small))

devtools::use_data(deseq_small, overwrite = TRUE, compress = "xz")
