\name{getIsoEIC.mz}
\alias{getIsoEIC.mz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get EIC of isotopic features.
}
\description{
Get EIC of isotopic features.
}
\usage{
getIsoEIC.mz(raw, fmz, nmax = 4, ppm = 50, rtrange = c(0, Inf),
                 charge = 1, C13 = 3, H2 = 0, O18 = 0, N15 = 0, S34 =
                 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
Result obtain by \code{LoadData}
}
  \item{fmz}{
The m/z of targeted compound. if the formula is given.
}
  \item{nmax}{
n, where at maximum M+n isotopes are count.
}
  \item{ppm}{
m/z range of EIC, given as ppm.
}
  \item{rtrange}{
Rt range of EIC.
}
  \item{charge}{
The number of charge of ion
}
  \item{C13}{
The maximum number of C13 atom.
}
  \item{H2}{
The maximum number of H2 atom.
}
  \item{O18}{
The maximum number of O18 atom.
}
  \item{N15}{
The maximum number of N15 atom.
}
  \item{S34}{
The maximum number of S34 atom.
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

eics <- getIsoEIC.mz(raw, fmz= 200.04)
}

