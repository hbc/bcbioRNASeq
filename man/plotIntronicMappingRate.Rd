% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R,
%   R/plotIntronicMappingRate-methods.R
\name{plotIntronicMappingRate}
\alias{plotIntronicMappingRate}
\alias{plotIntronicMappingRate,bcbioRNASeq-method}
\title{Plot intronic mapping rate}
\usage{
plotIntronicMappingRate(object, ...)

\S4method{plotIntronicMappingRate}{bcbioRNASeq}(
  object,
  interestingGroups = NULL,
  limit = 0.2,
  labels = list(title = "Intronic mapping rate", subtitle = NULL, sampleAxis = NULL,
    metricAxis = "intronic mapping rate (\%)"),
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
The majority of reads should map to exons and not introns.
}
\note{
Updated 2022-05-09.
}
\examples{
data(bcb)

## bcbioRNASeq ====
plotIntronicMappingRate(bcb)
}
\author{
Michael Steinbaugh, Rory Kirchner, Victor Barrera
}
