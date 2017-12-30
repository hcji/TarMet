\name{getEIC}
\alias{getEIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get EICs from raw data.
}
\description{
Get EICs from raw data.
}
\usage{
getEIC(raw, rtrange, mzrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
Result obtain by \code{LoadData} 
}
  \item{rtrange}{
Range of retention time
}
  \item{mzrange}{
Range of m/z
}
}

\examples{
# not run

# library(mtbls2)
# filepath <- file.path(find.package("mtbls2"), "mzData")
# files <- file.path(filepath, list.files(filepath, recursive = TRUE))[1:3]
# file1 <- files[1]
# raw <- LoadData(file1)

# eic <- getEIC(raw, rtrange = c(0, Inf), mzrange = c(201.03, 201.05))
# plot(eic$rt, eic$intensity, 'l')
}