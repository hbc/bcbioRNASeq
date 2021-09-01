.pkgName <- packageName()
.pkgVersion <- packageVersion(.pkgName)



.levels <- c("genes", "transcripts")

.assays <- "counts"
.tximportAssays <- c(.assays, "tpm")
.featureCountsAssays <- .assays

.tximportCallers <- c("salmon", "kallisto", "sailfish")
.featureCountsCallers <- c("star", "hisat2")
.callers <- c(.tximportCallers, .featureCountsCallers)

## DESeqTransform types.
.dt <- c("vst", "rlog")
.deseqAssays <- c("fpkm", .dt)

## "sf" denotes size-factor adjusted counts. This corresponds to
## `normalized = TRUE`, and is recommended by default for plots.
.normalized <- c("tpm", "sf", .deseqAssays, "tmm", "rle")

.legacyMetricsCols <- c("name", "x53Bias")



#' Cache URL
#'
#' @keywords internal
#' @export
#'
#' @examples
#' bcbioRNASeqTestsURL
bcbioRNASeqTestsURL <- paste0(
    "http://r.acidgenomics.com/testdata/bcbiornaseq/",
    "v", .pkgVersion$major, ".", .pkgVersion$minor  # nolint
)
