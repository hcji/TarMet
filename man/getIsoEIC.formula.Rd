\name{getIsoEIC.formula}
\alias{getIsoEIC.formula}

\title{
Get EIC of isotopic features.
}
\description{
Get EIC of isotopic features.
}
\usage{
getIsoEIC.formula(raw, formula, nmax = 3, adduct = "M+H", ppm = 50, rtrange = c(0, Inf), threshold = 0.01)
}

\arguments{
  \item{raw}{
Result obtain by \code{LoadData}
}
  \item{formula}{
A chemical formula
}
  \item{nmax}{
n, where at maximum M+n isotopes are count.
}
  \item{adduct}{
The species of adduct.
}
  \item{ppm}{
m/z range of EIC, given as ppm.
}
  \item{rtrange}{
Rt range of EIC.
}
  \item{threshold}{
The relative intensity threshold of isotope, from 0 to 1.
}
}

\value{
  A list of EICs of isotopic features.
}

\examples{
library(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- file.path(filepath, list.files(filepath, recursive = TRUE))[1:3]
file1 <- files[1]
raw <- LoadData(file1)

eics <- getIsoEIC.formula(raw, formula = 'C16H12O5', adduct = "M+H")
}