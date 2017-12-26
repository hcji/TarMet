\name{getIsoEIC}
\alias{getIsoEIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get EIC of isotopic features.
}
\description{
Get EIC of isotopic features.
}
\usage{
getIsoEIC(raw, formula, fmz, nmax = 3, adduct = "M+H", ppm = 50, rtrange = c(0, Inf), threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
Result obtain by \code{LoadData}
}
  \item{formula}{
A chemical formula or \code{NULL}
}
  \item{fmz}{
The m/z of targetted compound. if the formula is given, please set as -1.
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
Rt raneg of EIC.
}
  \item{threshold}{
The relative intensity threshold of isotope, from 0 to 1.
}
}
