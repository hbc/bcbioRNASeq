## Used for bcbio pipeline checks.
## Updated 2019-07-23.
.dataHasVariation <- function(dds) {
    !all(rowSums(assay(dds) == assay(dds)[, 1L]) == ncol(dds))
}



## Updated 2019-07-23.
.ddsMsg <- function() {
    message(paste0(
        "Generating DESeqDataSet with DESeq2 ",
        packageVersion("DESeq2"), "."
    ))
}



## Updated 2019-07-23.
`.new,DESeqDataSet` <-  # nolint
    function(se) {
        .ddsMsg()
        assert(is(se, "SummarizedExperiment"))

        ## Assert that counts are gene level.
        level <- metadata(se)[["level"]]
        assert(isString(level))
        if (level != "genes") {
            stop("Gene-level counts are required.")  # nocov
        }

        ## Subset the assays. Average transcript length matrix should only be
        ## included when raw counts are from tximport and not length scaled.
        if (
            metadata(se)[["caller"]] %in% tximportCallers &&
            metadata(se)[["countsFromAbundance"]] == "no"
        ) {
            message("Including average transcript length matrix.")
            assayNames <- c("counts", "avgTxLength")
        } else {
            assayNames <- "counts"
        }
        assays(se) <- assays(se)[assayNames]

        ## DESeq2 requires integer counts.
        counts <- counts(se)
        counts <- round(counts, digits = 0L)
        mode(counts) <- "integer"
        counts(se) <- counts

        ## Subset and update metadata.
        m <- metadata(se)
        keep <- c(
            "version",
            "interestingGroups",
            "uploadDir",
            "caller",
            "countsFromAbundance",
            "organism",
            "genomeBuild",
            "ensemblRelease",
            "runDate"
        )
        m <- m[intersect(names(m), keep)]
        m[["date"]] <- Sys.Date()
        m[["sessionInfo"]] <- session_info()
        metadata(se) <- m

        ## Generate the DESeqDataSet.
        ## Using an empty design formula.
        dds <- DESeqDataSet(se = se, design = ~ 1L)
        validObject(dds)
        dds
    }
