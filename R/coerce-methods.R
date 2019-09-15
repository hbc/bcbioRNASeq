#' Methods for coercing an object to a class
#'
#' Force an object to belong to a class.
#'
#' @name coerce
#' @author Michael Steinbaugh
#' @importFrom methods coerce
#' @exportMethod coerce
#' @note Updated 2019-08-07.
#'
#' @section bcbioRNASeq to DESeqDataSet:
#'
#' 1. Coerce to `RangedSummarizedExperiment`.
#' 2. Round raw counts to `integer matrix`.
#' 3. Subset [`colData()`][SummarizedExperiment::colData] to include only clean
#'    factor columns. See [`sampleData()`][basejump::sampleData] for details.
#' 4. Simplify [`metadata()`][S4Vectors::metadata] to include only relevant
#'    information and updates `sessionInfo`.
#'
#' Note that gene-level counts are required. Alternatively,
#' [`summarizeToGene()`][tximport::summarizeToGene] can be called to convert
#' transcript-level counts to gene-level. By default, we're using length-scaled
#' TPM, so a corresponding average transcript length matrix isn't necessary. The
#' average transcript length matrix is only necessary when raw counts matrix
#' isn't scaled during tximport call (see `countsFromAbundance` in
#' [`tximport()`][tximport::tximport] documentation).
#'
#' @section bcbioRNASeq to DESeqTransform:
#'
#' 1. Coerce to `DESeqDataSet`.
#' 2. Call [DESeq2::DESeq()].
#' 3. Call [DESeq2::varianceStabilizingTransformation()].
#'
#' @section bcbioRNASeq to DGEList:
#'
#' 1. Obtain per-observation scaling factors for length, adjusted to avoid
#'    changing the magnitude of the counts.
#' 2. Computing effective library sizes from scaled counts, to account for
#'    composition biases between samples.
#' 3. Combine effective library sizes with the length factors, and calculate
#'    offsets for a log-link GLM.
#' 4. Call [edgeR::DGEList()]
#' 4. Apply offset matrix using [edgeR::scaleOffset()].
#'
#' @seealso
#' - [tximport::tximport()].
#' - [DESeq2::DESeqDataSetFromTximport()].
#' - [edgeR::DGEList()].
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq to DESeqDataSet ====
#' x <- as(bcb, "DESeqDataSet")
#' names(S4Vectors::mcols(x))
#' class(x)
#' show(x)
#'
#' ## bcbioRNASeq to RangedSummarizedExperiment ====
#' x <- as(bcb, "RangedSummarizedExperiment")
#' slotNames(x)
#' show(x)
#'
#' ## bcbioRNASeq to SummarizedExperiment ====
#' ## Coerce to RangedSummarizedExperiment first.
#' x <- as(bcb, "RangedSummarizedExperiment")
#' x <- as(x, "SummarizedExperiment")
#' class(x)
#' slotNames(x)
#' show(x)
NULL



## Updated 2019-07-23.
`coerce,bcbioRNASeq,DESeqDataSet` <-  # nolint
    function(from) {
        validObject(from)
        rse <- as(from, "RangedSummarizedExperiment")
        `new,DESeqDataSet`(se = rse)
    }



#' @rdname coerce
#' @name coerce,bcbioRNASeq,DESeqDataSet-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqDataSet",
    def = `coerce,bcbioRNASeq,DESeqDataSet`
)



## Updated 2019-07-23.
`coerce,bcbioRNASeq,DESeqTransform` <-  # nolint
    function(from) {
        validObject(from)
        dds <- as(from, "DESeqDataSet")
        ## Expect warning about empty design formula.
        dds <- suppressWarnings(DESeq(dds))
        validObject(dds)
        message("Applying variance stabilizing transformation.")
        dt <- varianceStabilizingTransformation(dds)
        validObject(dt)
        dt
    }



#' @rdname coerce
#' @name coerce,bcbioRNASeq,DESeqTransform-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqTransform",
    def = `coerce,bcbioRNASeq,DESeqTransform`
)



## Note that we're following the tximport recommendations here.
## Updated 2019-07-23.
`coerce,bcbioRNASeq,DGEList` <-  # nolint
    function(from) {
        validObject(from)
        message(sprintf(
            "Generating DGEList with edgeR %s.",
            packageVersion("edgeR")
        ))
        ## Raw counts (i.e. txi$counts)
        cts <- counts(from)
        ## Average transcript length (i.e. txi$length)
        normMat <- assays(from)[["avgTxLength"]]
        ## Obtain per-observation scaling factors for length, adjusted to avoid
        ## changing the magnitude of the counts.
        normMat <- normMat / exp(rowMeans(log(normMat)))
        normCts <- cts / normMat
        ## Computing effective library sizes from scaled counts, to account for
        ## composition biases between samples.
        effLib <- calcNormFactors(normCts) * colSums(normCts)
        ## Combine effective library sizes with the length factors, and
        ## calculate offsets for a log-link GLM.
        normMat <- sweep(x = normMat, MARGIN = 2L, STATS = effLib, FUN = "*")
        normMat <- log(normMat)
        ## Creating a DGEList object for use in edgeR.
        to <- DGEList(cts)
        to <- scaleOffset(to, offset = normMat)
        ## Note that tximport guide recommends `filterByExpr()` step but we're
        ## intentionally skipping that step here.
        assert(identical(dimnames(from), dimnames(to)))
        validObject(to)
        to
    }



#' @rdname coerce
#' @name coerce,bcbioRNASeq,DGEList-method
setAs(
    from = "bcbioRNASeq",
    to = "DGEList",
    def = `coerce,bcbioRNASeq,DGEList`
)
