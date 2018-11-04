# Calculate a numeric vector to define the colors
# -1: downregulated
#  0: not significant
#  1: upregulated
.addIsDECol <- function(
    data,
    testCol = "padj",
    alpha,
    lfcCol = "log2FoldChange",
    lfcThreshold = 0L
) {
    # test: P value or S value
    test <- data[[testCol]]
    # lfc: log2 fold change cutoff
    lfc <- data[[lfcCol]]
    isDE <- mapply(
        test = test,
        lfc = lfc,
        FUN = function(test, lfc) {
            if (any(is.na(c(test, lfc)))) {
                # nonsignificant
                0L
            } else if (test < alpha & lfc > lfcThreshold) {
                # upregulated
                1L
            } else if (test < alpha & lfc < -lfcThreshold) {
                # downregulated
                -1L
            } else {
                0L
            }
        },
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
    )
    isDE <- as.factor(isDE)
    data[["isDE"]] <- isDE
    data
}



# Used for bcbio pipeline checks.
.dataHasVariation <- function(dds) {
    !all(rowSums(assay(dds) == assay(dds)[, 1L]) == ncol(dds))
}



.ddsMsg <- function() {
    message(paste0(
        "Generating DESeqDataSet with DESeq2 ",
        packageVersion("DESeq2"), "."
    ))
}



# Get differential expressed genes (DEGs) from DESeqResults table.
# Note that we're not sorting the identifiers here by LFC or P value.
# It's just performing a simple subset to get the identifiers as a character.
.deg <- function(
    object,
    alpha = NULL,
    lfcThreshold = NULL,
    direction = c("both", "up", "down")
) {
    assert_is_all_of(object, "DESeqResults")
    if (is.null(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    assertIsAlpha(alpha)
    if (is.null(lfcThreshold)) {
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
    }
    assert_is_a_number(lfcThreshold)
    assert_all_are_non_negative(lfcThreshold)
    direction <- match.arg(direction)

    # Define symbols to use in dplyr calls below.
    padj <- sym("padj")
    lfc <- sym("log2FoldChange")

    # Coerce to minimal tibble.
    data <- as(object, "tbl_df")
    cols <- c("rowname", "log2FoldChange", "padj")
    assert_is_subset(cols, colnames(data))
    data <- select(data, !!!syms(cols))

    # Apply alpha cutoff.
    data <- filter(data, !!padj < !!alpha)

    # Apply LFC threshold cutoff.
    if (lfcThreshold > 0L) {
        data <- filter(
            data,
            !!lfc > UQ(lfcThreshold) | !!lfc < -UQ(lfcThreshold)
        )
    }

    # Apply directional filtering.
    if (direction == "up") {
        data <- filter(data, !!lfc > 0L)
    } else if (direction == "down") {
        data <- filter(data, !!lfc < 0L)
    }

    deg <- pull(data, "rowname")
    message(paste(length(deg), "differentially expressed genes detected."))
    deg
}



.matchResults <- function(
    object,
    results,
    lfcShrink = FALSE
) {
    assert_is_all_of(object, "DESeqAnalysis")
    # Default to using the first contrast, for convenience.
    if (missing(results)) {
        results <- 1L
    }
    assert_is_scalar(results)
    assert_is_a_bool(lfcShrink)
    if (isTRUE(lfcShrink)) {
        slotName <- "lfcShrink"
    } else {
        slotName <- "results"
    }
    results <- slot(object, name = slotName)[[results]]
    assert_is_all_of(results, "DESeqResults")

    # Inform the user about which data we're using.
    msg <- paste("DESeqResults:", contrastName(results))
    if (isTRUE(lfcShrink)) {
        msg <- paste(msg, "(shrunken LFC)")
    }
    message(msg)

    results
}



.new.DESeqDataSet <-  # nolint
    function(se) {
        .ddsMsg()
        assert_that(is(se, "SummarizedExperiment"))

        # Assert that counts are gene level.
        level <- metadata(se)[["level"]]
        assert_is_a_string(level)
        if (level != "genes") {
            stop("Gene-level counts are required.")
        }

        # Subset the assays. Average transcript length matrix should only be
        # included when raw counts are from tximport and not length scaled.
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

        # DESeq2 requires integer counts.
        counts <- counts(se)
        counts <- round(counts, digits = 0L)
        mode(counts) <- "integer"
        counts(se) <- counts

        # Subset and update metadata.
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

        # Generate the DESeqDataSet.
        # Using an empty design formula.
        dds <- DESeqDataSet(se = se, design = ~ 1L)
        validObject(dds)
        dds
    }



.new.DESeqDataSetFromMatrix <-  # nolint
    function(countData) {
        .ddsMsg()
        assert_is_matrix(countData)
        # Integer counts are required.
        countData <- round(countData, digits = 0L)
        mode(countData) <- "integer"
        colData <- DataFrame(row.names = colnames(countData))
        dds <- DESeqDataSetFromMatrix(
            countData = countData,
            colData = colData,
            design = ~ 1L
        )
        validObject(dds)
        dds
    }



.transformCountsAxisLabel <- function(object) {
    paste(.transformType(object), "counts (log2)")
}



.transformType <- function(object) {
    assert_is_all_of(object, "DESeqTransform")
    if ("rlogIntercept" %in% colnames(mcols(object))) {
        "rlog"
    } else {
        # varianceStabilizingTransformation.
        "vst"
    }
}
