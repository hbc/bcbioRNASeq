% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotGenderMarkers-methods.R
\name{plotGenderMarkers}
\alias{plotGenderMarkers}
\alias{plotGenderMarkers,bcbioRNASeq-method}
\title{Plot sexually dimorphic gender marker genes}
\usage{
plotGenderMarkers(object, ...)

\S4method{plotGenderMarkers}{bcbioRNASeq}(
  object,
  normalized = c("tpm", "sf", "fpkm", "vst", "rlog", "tmm", "rle"),
  ...
)
}
\arguments{
\item{object}{Object.}

\item{normalized}{\code{character(1)} or \code{logical(1)}.
Normalization method to apply:
\itemize{
\item \code{FALSE}: Raw counts.
When using a \href{https://bioconductor.org/packages/tximport/}{tximport}-compatible caller, these are length scaled
by default (see \code{countsFromAbundance} argument).
When using a \href{https://doi.org/10.1093/bioinformatics/btt656}{featureCounts}-compatible caller, these are \code{integer}.
}

\href{https://bioconductor.org/packages/tximport/}{tximport} caller-specific normalizations:
\itemize{
\item \code{"tpm"}: \strong{T}ranscripts \strong{p}er \strong{m}illion.
}

Additional gene-level-specific normalizations:
\itemize{
\item \code{TRUE} / \code{"sf"}: \strong{S}ize \strong{f}actor (i.e. library size) normalized
counts.\cr
See \code{DESeq2::sizeFactors} for details.
\item \code{"fpkm"}: \strong{F}ragments \strong{p}er \strong{k}ilobase per \strong{m}illion mapped
fragments.\cr
Requires \code{fast = FALSE} in \code{bcbioRNASeq()} call and gene annotations in
\code{rowRanges()} with defined \code{width()}.\cr
See \code{DESeq2::fpkm()} for details.
\item \code{"vst"}: \strong{V}ariance-\strong{s}tabilizing \strong{t}ransformation (log2).\cr
Requires \code{fast = FALSE} to be set during \code{bcbioRNASeq()} call.\cr
See \code{DESeq2::varianceStabilizingTransformation()} for more information.
\item \code{"tmm"}: \strong{T}rimmed \strong{m}ean of \strong{M}-values.\cr
Calculated on the fly.\cr
See \code{edgeR::calcNormFactors()} for details.
\item \code{"rle"}: \strong{R}elative \strong{l}og \strong{e}xpression transformation.\cr
Calculated on the fly.\cr
See \code{relativeLogExpression()} for details.
\item \code{"rlog"}: \emph{Deprecated}.
\strong{R}egularized \strong{log} transformation (log2).\cr
No longer calculated automatically during \code{bcbioRNASeq()} call, but may
be defined in legacy objects.\cr
See \code{DESeq2::rlog()} for details.\cr
Note that VST is more performant and now recommended by default instead.
}

Note that \code{logical(1)} support only applies to \code{counts()}. Other
functions in the package require \code{character(1)} and use \code{match.arg()}
internally.}

\item{...}{Passthrough to \code{SummarizedExperiment} method defined in AcidPlots.
See \code{AcidPlots::plotGenderMarkers()} for details.}
}
\value{
\code{ggplot}.
}
\description{
This is a convenience function that quickly plots known sexually dimorphic
genes, and overall gene experssion for the X and Y chromosomes.
}
\note{
Updated 2022-03-07.
}
\examples{
data(bcb)

## bcbioRNASeq ====
## Simulate expression of sexually dimorphic gender marker genes.
rownames(bcb)[seq_len(2L)] <- c("ENSMUSG00000086503", "ENSMUSG00000069045")
plotGenderMarkers(bcb)
}
\seealso{
\code{plotCounts()}.
}
\author{
Michael Steinbaugh
}
