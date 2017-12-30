\name{PAFFT}
\alias{PAFFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alignment based on PAFFT.
}
\description{
Alignment based on PAFFT.
}
\usage{
PAFFT(spectra, reference, segSize, shift)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spectra}{
The spectras to be aligned.
}
  \item{reference}{
The reference.
}
  \item{segSize}{
The size of segment.
}
  \item{shift}{
Maxmimum shift.
}
}

\examples{
data(exampleMS, package = 'MassSpecWavelet')
spectra = t(c(exampleMS, rep(0,50)))
reference = t(c(rep(0,50), exampleMS))
res = PAFFT(spectra, reference, 1000, 100)
}