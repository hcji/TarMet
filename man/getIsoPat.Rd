\name{getIsoPat}
\alias{getIsoPat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the theoretical isotope pattern.
}
\description{
Get the theoretical isotope pattern.
}
\usage{
getIsoPat(formula, d, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A chemical formula.
}
  \item{d}{
The species of adduct
}
  \item{threshold}{
The relative intensity threshold, from 0 to 1.
}
}

\examples{
pattern <- getIsoPat('C16H12O5', 'M+H', 0.2)
}