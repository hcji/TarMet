\name{getArea}
\alias{getArea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get peak area of a peak.
}
\description{
 Get peak area of a peak.
}
\usage{
getArea(eics, rtmin, rtmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eics}{
EICs obtained by \code{getIsoEIC} function.
}
  \item{rtmin}{
Start of retention time.
}
  \item{rtmax}{
End of retention time.
}
}

\value{
  the peak areas.
}

\examples{
library(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- file.path(filepath, list.files(filepath, recursive = TRUE))[1:3]
file1 <- files[1]
raw <- LoadData(file1)
eics <- getIsoEIC(raw, formula = NULL, fmz= 200.04, adduct = "M+H")

areas <- getArea(eics, rtmin = 0, rtmax = Inf)
}