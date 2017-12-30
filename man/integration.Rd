\name{integration}
\alias{integration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{To compute the peak area under the each fitting peaks by trapezoidal integration.}
\description{
 To compute the peak area under the each fitting peaks by trapezoidal integration.
}
\usage{
integration(x,yf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The values of x-axis}
  \item{yf}{The yf is a collection of each fitting peaks.}
}

\value{
  the peak area obtained by integration.
}

\examples{
data(exampleMS, package = 'MassSpecWavelet')
integration(1:length(exampleMS), exampleMS)
}