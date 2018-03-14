library(pryr)
library(devtools)
library(DESeq2)
load_all()

dds_small <- as(bcb_small, "DESeqDataSet")
design(dds_small) <- ~ treatment
dds_small <- DESeq(dds_small)
validObject(dds_small)
stopifnot(design(dds_small) == ~ treatment)
stopifnot(identical(
    resultsNames(dds_small),
    c("Intercept", "treatment_folic_acid_vs_control")
))

rld_small <- rlog(dds_small)
validObject(rld_small)

res_small <- results(
    dds_small,
    contrast = c(
        factor = "treatment",
        numerator = "folic_acid",
        denominator = "control"
    )
)
validObject(res_small)

use_data(dds_small, rld_small, res_small, overwrite = TRUE, compress = "xz")
