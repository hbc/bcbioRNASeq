% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/deprecated.R
\name{plotDegHeatmap}
\alias{plotDegHeatmap}
\alias{deprecated}
\alias{camel}
\alias{plotCountsPerGene}
\alias{plotDEGHeatmap}
\alias{plotDegHeatmap,missing-method}
\alias{plotGeneSaturation}
\alias{plotGenesDetected}
\alias{plotMA}
\alias{plotMeanAverage}
\alias{plotMeanSD}
\alias{plotPCA}
\alias{plotPCACovariates}
\alias{plotQC}
\alias{plotRRNAMappingRate}
\alias{prepareRNASeqTemplate}
\alias{topTables}
\alias{writeCounts}
\title{Deprecated functions}
\usage{
plotDegHeatmap(object, ...)

camel(...)

plotCountsPerGene(object, ...)

plotDEGHeatmap(...)

\S4method{plotDegHeatmap}{missing}(object, results, counts, ...)

plotGeneSaturation(object, ...)

plotGenesDetected(object, ...)

plotMA(...)

plotMeanAverage(...)

plotMeanSD(...)

plotPCA(...)

plotPCACovariates(...)

plotQC(...)

plotRRNAMappingRate(...)

prepareRNASeqTemplate(...)

topTables(...)

writeCounts(..., dir = getOption(x = "acid.export.dir", default = getwd()))
}
\arguments{
\item{...}{Additional arguments.}

\item{object}{Object.}

\item{counts}{\code{matrix}.
Count matrix. Normalized counts are recommended.}

\item{dir}{\code{character(1)}.
Directory path.}

\item{`results`}{\code{DESeqResults.}}

\item{`counts`}{\code{DESeqTransform}.}
}
\value{
\code{.Deprecated}.
}
\description{
Deprecated functions
}
\examples{
a <- function(...) {
    .Deprecated("b")
    b(...)
}

b <- function(x) {
    x + 1L
}

withCallingHandlers(
    expr = a(1L),
    warning = function(w) {
        print(w)
        invokeRestart("muffleWarning")
    }
)
}
\seealso{
\link[base:Deprecated]{Deprecated}.
}
\keyword{internal}
