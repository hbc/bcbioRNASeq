library(devtools)
library(tidyverse)
load_all()
meta <- colData(bcb_small) %>%
    mutate(
        # Add required `fileName` column
        fileName = paste0(sampleID, ".fastq.gz"),
        # Ensure `sampleID` and `sampleName` aren't defined
        sampleID = NULL,
        sampleName = NULL
    ) %>%
    select(fileName, everything())
write_csv(meta, path = "tests/testthat/sample_metadata.csv")
