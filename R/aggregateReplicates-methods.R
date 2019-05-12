#' @name aggregateReplicates
#' @inherit bioverbs::aggregateReplicates
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `RangedSummarizedExperiment`.
#'
#' @examples
#' bcb <- bcb_small
#' colnames(bcb)
#' # Assign groupings into `aggregate` column of `colData()`
#' aggregate <- as.factor(sub("^([a-z]+)_.*", "\\1", colnames(bcb)))
#' names(aggregate) <- colnames(bcb)
#' aggregate
#' bcb$aggregate <- aggregate
#' aggregateReplicates(bcb)
NULL



#' @rdname aggregateReplicates
#' @name aggregateReplicates
#' @importFrom bioverbs aggregateReplicates
#' @usage aggregateReplicates(object, ...)
#' @export
NULL



aggregateReplicates.bcbioRNASeq <-  # nolint
    function(object) {
        validObject(object)

        metadata <- metadata(object)
        colData <- colData(object)
        assert_is_subset("aggregate", colnames(colData))
        assert_is_factor(colData[["aggregate"]])

        # This step will replace the `sampleName` column with the `aggregate`
        # column metadata.
        remap <- colData %>%
            as.data.frame() %>%
            rownames_to_column("sampleID") %>%
            select(!!!syms(c("sampleID", "aggregate"))) %>%
            mutate(sampleIDAggregate = makeNames(
                !!sym("aggregate"), unique = FALSE
            )) %>%
            select(-!!sym("aggregate")) %>%
            arrange(!!!syms(c("sampleID", "sampleIDAggregate"))) %>%
            mutate_all(as.factor) %>%
            mutate_all(droplevels)

        # Message the new sample IDs
        newIDs <- unique(remap[["sampleIDAggregate"]])
        message(paste("New sample IDs:", toString(newIDs)))

        groupings <- factor(remap[["sampleIDAggregate"]])
        names(groupings) <- remap[["sampleID"]]

        # Assays ---------------------------------------------------------------
        message("Aggregating counts")
        counts <- aggregateReplicates(counts(object), groupings = groupings)
        assert_are_identical(sum(counts), sum(counts(object)))

        # Column data ----------------------------------------------------------
        # Return minimal metadata with `sampleName` column only
        expected <- length(levels(colData[["aggregate"]]))
        colData <- colData %>%
            as.data.frame() %>%
            mutate(sampleName = !!sym("aggregate")) %>%
            select(!!sym("sampleName")) %>%
            mutate_all(as.factor) %>%
            unique() %>%
            as("DataFrame")
        if (!identical(nrow(colData), expected)) {
            stop("Failed to aggregate sample metadata uniquely")  # nocov
        }
        rownames(colData) <- makeNames(colData[["sampleName"]])
        assert_are_identical(colnames(counts), rownames(colData))

        # Return ---------------------------------------------------------------
        SummarizedExperiment(
            assays = list(counts = counts),
            colData = colData,
            rowRanges = rowRanges(object)
        )
    }



#' @rdname aggregateReplicates
#' @export
setMethod(
    f = "aggregateReplicates",
    signature = signature("bcbioRNASeq"),
    definition = aggregateReplicates.bcbioRNASeq
)
