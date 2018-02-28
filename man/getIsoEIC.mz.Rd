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
getIsoEIC.mz(raw, fmz, nmax = 4, ppm = 50, rtrange = c(0, Inf), charge = 1, elements = NULL, isotope = NULL, numbers = NULL)
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

  \item{elements}{
The type of elements
}
  \item{isotope}{
The isotopes of elements to be considered
}
  \item{numbers}{
n, where at maximum M+n isotopes are count.
}
}
\value{
  A list of EICs of isotopic features.
}