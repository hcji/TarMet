\name{getIsoPeaks}
\alias{getIsoPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get peaks of isotopic EICs.
}
\description{
Get peaks of isotopic EICs.
}
\usage{
getIsoPeaks(eics, SNR.Th = 4, peakScaleRange = 5, peakThr = 0, userDefTime = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eics}{
EICs obtained by \code{getIsoEIC} function.
}
  \item{SNR.Th}{
Threshold of SNR.
}
  \item{peakScaleRange}{
The minimum CWT scale range of peaks.
}
  \item{peakThr}{
The minimum intensity of peaks
}
  \item{userDefTime}{
User defined retention time ranges of a select peak.
}
}

\value{
  the peak information of isotopic peaks.
}

\examples{
library(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- file.path(filepath, list.files(filepath, recursive = TRUE))[1:3]

file1 <- files[1]
raw <- LoadData(file1)
eics <- getIsoEIC.mz(raw, fmz= 200.04) 

peaks <- getIsoPeaks(eics, SNR.Th = 4, peakScaleRange = 5)
}
