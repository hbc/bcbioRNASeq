% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R,
%   R/relativeLogExpression-methods.R
\name{relativeLogExpression}
\alias{relativeLogExpression}
\alias{relativeLogExpression,matrix-method}
\alias{relativeLogExpression,SummarizedExperiment-method}
\title{Relative log expression}
\usage{
relativeLogExpression(object, ...)

\S4method{relativeLogExpression}{matrix}(object)

\S4method{relativeLogExpression}{SummarizedExperiment}(object)
}
\arguments{
\item{object}{Object.}

\item{...}{Additional arguments.}
}
\value{
\code{matrix}.
}
\description{
Relative log expression
}
\note{
Updated 2020-01-20.
}
\examples{
## bcbioRNASeq ====
data(bcb)
relativeLogExpression(bcb)
}
\references{
Anders and Huber (2010).
}
\seealso{
\code{edgeR::calcNormFactors()}.
}
\author{
Lorena Pantano, Michael Steinbaugh
}
