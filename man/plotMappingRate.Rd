% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotMappingRate-methods.R
\name{plotMappingRate}
\alias{plotMappingRate}
\alias{plotMappingRate,bcbioRNASeq-method}
\title{Plot mapping rate}
\usage{
plotMappingRate(object, ...)

\S4method{plotMappingRate}{bcbioRNASeq}(
  object,
  interestingGroups = NULL,
  limit = 0.7,
  labels = list(title = "Mapping rate", subtitle = NULL, sampleAxis = NULL, metricAxis =
    "mapping rate (\%)"),
  flip = getOption(x = "acid.flip", default = TRUE)
)
}
\arguments{
\item{object}{Object.}

\item{interestingGroups}{\code{character}.
Groups of interest to use for visualization.
Corresponds to factors describing the columns of the object.}

\item{limit}{\code{numeric(1)}.
Threshold limit.}

\item{labels}{\code{list}.
ggplot2 labels.
See \code{ggplot2::labs()} for details.}

\item{flip}{\code{logical(1)}.
Flip x and y axes. Recommended for plots containing many samples.}

\item{...}{Additional arguments.}
}
\value{
\code{ggplot}.
}
\description{
The genomic mapping rate represents the percentage of reads mapping to the
reference genome. Low mapping rates are indicative of sample contamination,
poor sequencing quality or other artifacts.
}
\note{
Updated 2022-05-09.
}
\examples{
data(bcb)

## bcbioRNASeq ====
plotMappingRate(bcb)
}
\author{
Michael Steinbaugh, Rory Kirchner, Victor Barrera
}
