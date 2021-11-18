LoadMSMS <- function(filename)
{
  # filename <- "D:/data/IROA_PLATE_pos_3_C.mzML"
  splitname <- strsplit(filename,"\\.")[[1]]
  if(tolower(splitname[length(splitname)]) == "cdf"){
    msobj <- openMSfile(filename,backend="netCDF")
  }else if (tolower(splitname[length(splitname)]) == "mzml"){
    msobj <- openMSfile(filename,backend="pwiz")
  }else{
    msobj <- openMSfile(filename,backend="Ramp")
  }
  
  peakInfo <- peaks(msobj)
  headerInfo <- header(msobj)
  whMS2 <- which(headerInfo$msLevel==2)
  
  if (length(whMS2) == 0){
    return(NULL)
  }
  
  peakInfo <- peakInfo[whMS2]
  peakInfo <- lapply(peakInfo, function(spectrum) {
    keep <- spectrum[,2] > 1e-6
    output <- as.data.frame(spectrum[keep,,drop = FALSE])
    colnames(output) <- c('mz','intensity')
    return(output)
  })
  
  scanTime <- round(headerInfo$retentionTime[whMS2],3)
  precursorMz <- round(headerInfo$precursorMZ[whMS2],3)
  precursorIntensity <- round(headerInfo$precursorIntensity[whMS2])
  headerInfo <- data.frame(scanTime=scanTime, precursorMz=precursorMz, precursorIntensity=precursorIntensity)
  
  return(list(header=headerInfo, peaks = peakInfo))
}


getMSMS <- function(header, peaks, precursorMzRange, precursorRtRange, ppm = 100){
  # precursorRtRange <- c(185.16, 367.16)
  # precursorMzRange <- c(162.1092, 162.1157)
  whMS2 <- which(header$precursorMz >= precursorMzRange[1] &
                 header$precursorMz <= precursorMzRange[2] &
                 header$scanTime >= precursorRtRange[1] &
                 header$scanTime <= precursorRtRange[2])
  whPeaks <- peaks[whMS2]
  spectrums <- list()
  if(length(whPeaks) > 0){
    for(p in whPeaks){
      s <- new("Spectrum2", rt = 1, precursorMz = mean(precursorMzRange),
                mz = p[,1], intensity = p[,2])
      spectrums <- c(spectrums, s)
    }
    ms2 <- MSnbase::consensusSpectrum(spectrums, minProp = 0.5, ppm = ppm)
    ms2 <- cbind(MSnbase::mz(ms2), cbind(MSnbase::intensity(ms2)))
    colnames(ms2) <- c('mz', 'intensity')
    return(ms2)
  } else {
    return(NULL)
  }
}

getMSVec <- function(ms, dmz=0.01){
  idx <- seq(10, 2000, dmz)
  vec <- rep(0, length(idx))
  for (i in 1:nrow(ms)){
    m <- ms[i,1]
    w <- which.min(abs(idx - m))
    vec[w] <- vec[w] + ms[i,2]
  }
  vec <- vec / max(vec)
  return(vec)
}

compareMS <- function(ref, que, dmz){
  refv <- getMSVec(ref, dmz)
  quev <- getMSVec(que, dmz)
  score <- refv %*% quev / sqrt((refv %*% refv) * (quev %*% quev))
  score <- as.numeric(score)
  return(score)
}
