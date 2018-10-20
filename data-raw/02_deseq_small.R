#' Example DESeq2 differential expression analysis
#' 2018-10-20

# Restrict to 2 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
library(pryr)
limit <- structure(2e6, class = "object_size")

library(DESeq2)

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
res_shrunken <- lfcShrink(
    dds = dds,
    coef = 2L,
    res = res
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
    FUN = object_size,
    FUN.VALUE = numeric(1L)
)
object_size(deseq_small)
stopifnot(object_size(bcb) < limit)

# Check that object is valid.
stopifnot(is(deseq_small, "DESeqAnalysis"))
stopifnot(validObject(deseq_small))

devtools::use_data(deseq_small, overwrite = TRUE, compress = "xz")
