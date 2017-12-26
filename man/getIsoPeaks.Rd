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