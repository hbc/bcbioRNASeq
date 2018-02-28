devtools::load_all()
library(tidyverse)

loadRemoteData("http://bcbiornaseq.seq.cloud/f1000v2/bcb_all.rda")

# Select only the control and day 7 folic acid samples
bcb <- selectSamples(bcb_all, day = c(0L, 7L))

# Fix metadata slots that take up too much disk space
metadata(bcb)$annotable <- metadata(bcb)$annotable[, c("ensgene", "symbol")]
metadata(bcb)$bcbioCommandsLog <- ""
metadata(bcb)$bcbioLog <- ""
metadata(bcb)$tx2gene <- head(metadata(bcb)$tx2gene, 2L)

# Select only the top 500 most abundant genes
tpm <- tpm(bcb)
genes <- rowSums(tpm) %>%
    sort(decreasing = TRUE) %>%
    head(100L) %>%
    names()
# Make sure dimorphic gender markers are included
dimorphic <- genderMarkers$musMusculus %>%
    filter(include == TRUE) %>%
    pull(ensgene) %>%
    sort() %>%
    .[. %in% rownames(bcb)]
genes <- c(genes, dimorphic) %>%
    unique() %>%
    sort() %>%
    head(100L)
bcb <- bcb[genes, ]

dds <- bcbio(bcb, "DESeqDataSet")
design(dds) <- formula(~day)
dds <- DESeq(dds)
rld <- rlog(dds)
res <- results(
    dds,
    contrast = c(
        factor = "day",
        numerator = 7L,
        denominator = 0L
    )
)

saveData(
    bcb, dds, rld, res,
    dir = file.path("inst", "extdata"),
    compress = "xz",
    overwrite = TRUE)
