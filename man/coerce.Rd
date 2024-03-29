% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/coerce-methods.R
\name{as.DESeqDataSet}
\alias{as.DESeqDataSet}
\alias{as.DESeqTransform}
\alias{as.DGEList}
\alias{coerce}
\alias{as.DESeqDataSet,bcbioRNASeq-method}
\alias{as.DESeqTransform,bcbioRNASeq-method}
\alias{as.DGEList,bcbioRNASeq-method}
\alias{coerce,bcbioRNASeq,DESeqDataSet-method}
\alias{coerce,bcbioRNASeq,DESeqTransform-method}
\alias{coerce,bcbioRNASeq,DGEList-method}
\title{Methods for coercing an object to a class}
\usage{
as.DESeqDataSet(x, ...)

as.DESeqTransform(x, ...)

as.DGEList(x, ...)

\S4method{as.DESeqDataSet}{bcbioRNASeq}(x, quiet = FALSE)

\S4method{as.DESeqTransform}{bcbioRNASeq}(x, quiet = FALSE)

\S4method{as.DGEList}{bcbioRNASeq}(x, quiet = FALSE)
}
\arguments{
\item{x}{Object.}

\item{quiet}{\code{logical(1)}.
Perform command quietly, suppressing messages.}

\item{...}{Additional arguments.}
}
\value{
Modified object, of desired coercion type.
}
\description{
Force an object to belong to a class.
}
\note{
Updated 2022-03-07.
}
\section{bcbioRNASeq to DESeqDataSet}{

\enumerate{
\item Coerce to \code{RangedSummarizedExperiment}.
\item Round raw counts to \verb{integer matrix}.
\item Subset \code{colData()} to include only clean
factor columns. See \code{sampleData()} for details.
\item Simplify \code{metadata()} to include only relevant information and
updates \code{sessionInfo}.
}

Note that gene-level counts are required. Alternatively,
\code{tximport::summarizeToGene()} can be called to convert transcript-level
counts to gene-level. By default, we're using length-scaled TPM, so a
corresponding average transcript length matrix isn't necessary. The average
transcript length matrix is only necessary when raw counts matrix isn't
scaled during tximport call (see \code{countsFromAbundance} in
\code{tximport::tximport()} documentation).
}

\section{bcbioRNASeq to DESeqTransform}{

\enumerate{
\item Coerce to \code{DESeqDataSet}.
\item Call \code{DESeq2::DESeq()}.
\item Call \code{DESeq2::varianceStabilizingTransformation()}.
}
}

\section{bcbioRNASeq to DGEList}{


When \code{countsFromAbundance = "lengthScaledTPM"} (default):
\enumerate{
\item Call \code{edgeR::DGEList()}.
}

When \code{countsFromAbundance = "no"}:
\enumerate{
\item Call \code{edgeR::DGEList()}.
\item Obtain per-observation scaling factors for length, adjusted to avoid
changing the magnitude of the counts.
\item Computing effective library sizes from scaled counts, to account for
composition biases between samples.
\item Combine effective library sizes with the length factors, and calculate
offsets for a log-link GLM.
\item Apply offset matrix using \code{edgeR::scaleOffset()}.
}
}

\examples{
data(bcb)

## bcbioRNASeq to DESeqDataSet ====
dds <- as.DESeqDataSet(bcb)
class(dds)

## bcbioRNASeq to DESeqTransform ====
dt <- as.DESeqTransform(bcb)
class(dt)
}
\seealso{
\itemize{
\item \code{tximport::tximport()}.
\item \code{DESeq2::DESeqDataSetFromTximport()}.
\item \code{edgeR::DGEList()}.
}
}
\author{
Michael Steinbaugh
}
