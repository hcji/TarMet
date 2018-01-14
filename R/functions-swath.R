LoadSwath <- function(filename) {
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
  
  whMS2 <- headerInfo$msLevel==2
  preMz_diff <- diff(headerInfo$precursorMZ[whMS2])
  preMz_win <- mean(preMz_diff[preMz_diff>0])+1
  
  peakInfo <- peakInfo[whMS2]
  peakInfo <- lapply(peakInfo, function(spectrum) {
    keep <- spectrum[,2] > 1e-6
    output <- as.data.frame(spectrum[keep,,drop = FALSE])
    colnames(output) <- c('mz','intensity')
    return(output)
  })
  
  scanTime <- round(headerInfo$retentionTime[whMS2],3)
  precursorMZ <- headerInfo$precursorMZ[whMS2]
  # close(msobj)
  
  return(list(path=filename, preMz_win=preMz_win, times=scanTime, precursorMZ = precursorMZ, peaks = peakInfo))
}

getSwathData <- function(swath_raw, tarMz){
  whPre <- abs(swath_raw$precursorMZ - tarMz) < 0.5*swath_raw$preMz_win
  peaks <- swath_raw$peaks[whPre]
  times <- swath_raw$times[whPre]
  
  return(list(times=times, peaks=peaks))
}

getSwathMz <- function(raw_swath, peak_position, ppm) {
  peak_position <- as.numeric(as.character(peak_position))
  scan <- which.min(abs(peak_position - (raw_swath$times)))
  SwathMz <- raw_swath$peaks[[scan]]$mz
  SwathInt <- raw_swath$peaks[[scan]]$intensity
  dels <- which(diff(SwathMz)/SwathMz[seq_along(SwathMz)[-1]]*10^6 < 0.5*ppm)
  if (length(dels)>0){
    SwathMz <- SwathMz[-dels]
  }
  return(SwathMz)
}

getSwathEICs <- function(raw_swath, swath_mz, ppm, rtrange = c(0, Inf)){
  mzranges <- do.call(rbind, lapply(swath_mz, function(mz){
    c(mz * (1 - ppm/2/10^6), mz * (1 + ppm/2/10^6))
  }))
  eics <- apply(mzranges, 1, function(mzrange){
    getEIC(raw_swath, rtrange, mzrange)
  })
  return (list(mzs=mzranges, eics=eics))
}

getPeakCorr <- function(eics_ms1, eics_ms2, tar_left, tar_right, ref_ind, align_seg, align_shift){
  eics_ms1 <- lapply(eics_ms1, function(f){
    f$eics[[1]]
  })
  
  width <- tar_right - tar_left
  rt_ref <- eics_ms1[[1]]$rt
  
  align_seg <- round(align_seg/mean(diff(rt_ref)))
  align_shift <- round(align_shift/mean(diff(rt_ref)))
  
  intensity_ms1 <- lapply(eics_ms1, function(f){
    new_intensity <- approx(f$rt, f$intensity, rt_ref, rule = 2)$y
    WhittakerSmooth(new_intensity, rep(1, length(new_intensity)), lambda = 3)
  })
  intensity_ms1 <- do.call(rbind, intensity_ms1)
  
  mz_frag <- eics_ms2[[1]]$mzs
  intensity_ms2 <- lapply(1:nrow(mz_frag), function(m){
    this <- lapply(eics_ms2, function(f){
      this.eic <- f$eics[[m]]
      this.intensity <- approx(this.eic$rt, this.eic$intensity, rt_ref, rule = 2)$y
      WhittakerSmooth(this.intensity, rep(1, length(this.intensity)), lambda = 3)
    })
    do.call(rbind, this)
  })
  
  intensity_ms1 <- PAFFT(intensity_ms1, intensity_ms1[ref_ind,], align_seg, align_shift)$alignedSpectrum
  intensity_ms2 <- lapply(seq_along(intensity_ms2), function(s){
    reference <- intensity_ms2[[s]][ref_ind,]
    spectra <- intensity_ms2[[s]]
    PAFFT(spectra, reference, align_seg, align_shift)$alignedSpectrum
  })
  
  tar_inds <- which(rt_ref>tar_left-0.5*width & rt_ref<tar_right+0.5*width)
  pearson_corr <- sapply(intensity_ms2, function(m){
    calPearsonCorr(intensity_ms1[, tar_inds], m[, tar_inds])
  })
  
  return(pearson_corr)
}

calPearsonCorr <- function(A,B){
  a <- A-mean(A)
  b <- B-mean(B)
  return(sum(a*b)/sqrt(sum(a*a)*sum(b*b)))
}
