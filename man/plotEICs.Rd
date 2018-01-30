\name{plotEICs}
\alias{plotEICs}

\title{
Plot EICs of objects.
}
\description{
Plot EICs of objects.
}
\usage{
plotEICs(eics, peaks = NULL, Names = NULL, rt=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eics}{
A list of eics.
}
  \item{peaks}{
Peak detection results or NULL
}
  \item{Names}{
Names of traces or NULL
}
  \item{rt}{
Range of retention time to be plot.
}
}

\value{
A EICs plot.
}

\examples{
library(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- file.path(filepath, list.files(filepath, recursive = TRUE))[1:3]
file1 <- files[1]
raw <- LoadData(file1)

eics <- getIsoEIC.formula(raw, formula = 'C16H12O5', adduct = "M+H")
plotEICs(eics)
}
