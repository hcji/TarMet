LoadData <- function(filename)
{
  library(mzR)
  splitname <- strsplit(filename,"\\.")[[1]]
  if(tolower(splitname[length(splitname)]) == "cdf")
  {
    msobj <- openMSfile(filename,backend="netCDF")
  }else{
    msobj <- openMSfile(filename,backend="Ramp")
  }
  
  peakInfo <- peaks(msobj)
  headerInfo <- header(msobj)
  whMS1 <- which(headerInfo$msLevel==1)
  peakInfo <- peakInfo[whMS1]
  
  peakInfo <- lapply(peakInfo, function(spectrum) {
    keep <- spectrum[,2] > 1e-6
    output <- as.data.frame(spectrum[keep,,drop = FALSE])
    colnames(output) <- c('mz','intensity')
    return(output)
  })
  
  scanTime <- round(headerInfo$retentionTime[whMS1],3)
  # close(msobj)
  
  return(list(path=filename, times=scanTime, peaks = peakInfo))
}

getEIC <- function(raw, rtrange, mzrange){
  scanrange <- findInterval(rtrange, raw$times)
  rt <- raw$times[(scanrange[1]+1) : scanrange[2]]
  intensity <- unlist(lapply((scanrange[1]+1) : scanrange[2], function(s){
    ind <- findInterval(mzrange, raw$peaks[[s]]$mz)
    if (ind[2] > ind[1]){
      return(sum(raw$peaks[[s]]$intensity[(ind[1]+1) : ind[2]]))
    } else {
      return(0)
    }
  }))
  return(list(rt = rt, intensity = intensity))
}

getIsoPat <- function(formula, d, threshold){
  data("isotopes", package = "enviPat")
  data("adducts", package = "enviPat")
  checked <- check_chemform(isotopes, formula)
  if (adduct > 0){
    if (adducts[d, 7] != 'FALSE'){
      chemform <- mergeform(checked[,2], adducts[d, 7])
    } else {
      chemform <- checked[, 2]
    }
    if (adducts[d, 8] != 'FALSE'){
      chemform <- subform(chemform, adducts[d, 8])
    }
  }
  checked <- check_chemform(isotopes, chemform)
  pattern <- isopattern(isotopes, checked$new_formula, threshold = threshold)[[1]][, c(1,2)]
  return(as.data.frame(pattern))
}

getIsoEIC <- function(raw, formula, fmz, adduct = 1, ppm = 50, rtrange = c(0, Inf), threshold = 0.01){
  
  if (fmz < 0){
    pattern <- getIsoPat(formula, adduct, threshold)
    nmax <- round(max(pattern[,1] - pattern[1,1]))
    mzs <- sapply(0:nmax, function(n){
      pai <- pattern[round(pattern[,1]- pattern[1,1])==n,]
      pai[which.max(pai[,2]) ,1]
    })
  } else {
    data("adducts", package = "enviPat")
    mzs <- (fmz + 0:2 * 1.0033) / abs(adducts[adduct, 3])
  }
  
  mzranges <- do.call(rbind, lapply(mzs, function(mz){
    mz <- mz + adducts$Mass[adduct]
    c(mz * (1 - ppm/2/10^6), mz * (1 + ppm/2/10^6))
  }))
  
  eics <- apply(mzranges, 1, function(mzrange){
    getEIC(raw, rtrange, mzrange)
  })
  
  return(list(mzs = mzranges, eics = eics))
}

getIsoPeaks <- function(eics, SNR.Th = 4, peakScaleRange = 5, peakThr = 0, userDefTime = NULL){
  eic <- eics$eics[[1]]
  MajorPeaks <- peakDetectionCWT(eic$intensity, SNR.Th = SNR.Th, peakScaleRange = peakScaleRange, peakThr = peakThr)
  PeakWidths <- widthEstimationCWT(eic$intensity, MajorPeaks$majorPeakInfo)
  
  Position <- MajorPeaks$majorPeakInfo$peakIndex
  Start <- sapply(2:length(PeakWidths), function(s){
    PeakWidths[[s]][1]
  })
  End <- sapply(2:length(PeakWidths), function(s){
    PeakWidths[[s]][length(PeakWidths[[s]])]
  })
  
  Area <- do.call(rbind, lapply(eics$eics, function(eic){
    sapply(1:length(Position), function(p){
      integration(eic$rt[Start[p]:End[p]], eic$intensity[Start[p]:End[p]] )
    })
  }))
  PeakArea <- as.data.frame(cbind(eics$mzs, Area))
  colnames(PeakArea) <- c('MzMin', 'MzMax', paste('Peak', 1:length(Position)))
  PeakInfo <- as.data.frame(cbind(paste('Peak', 1:length(Position)), eic$rt[Position], eic$rt[Start], eic$rt[End]))
  colnames(PeakInfo) <- c(' ', 'Position', 'Start', 'End')
  Index <- as.data.frame(cbind(Position, Start, End))
  
  return(list(PeakInfo = PeakInfo, PeakArea = PeakArea, Index = Index))
}

getArea <- function(eics, rtmin, rtmax){
  areas <- sapply(eics$eics, function(eic){
    Start <- findInterval(rtmin, eic$rt)+1
    End <- findInterval(rtmax, eic$rt)
    Area <- integration(eic$rt[Start:End], eic$intensity[Start:End])
    Area
  })
  
  return(areas)
}

airPLS <- function(x,lambda=10,differences=1, itermax=20){
  
  x = as.vector(x)
  m = length(x)
  w = rep(1,m)
  control = 1
  i = 1
  while(control==1){
    z = WhittakerSmooth(x,w,lambda,differences)
    d = x-z
    sum_smaller = abs(sum(d[d<0])) 
    if(sum_smaller<0.001*sum(abs(x))||i==itermax)
    {
      control = 0
    }
    w[d>=0] = 0
    w[d<0] = exp(i*abs(d[d<0])/sum_smaller)
    w[1] = exp(i*max(d[d<0])/sum_smaller)
    w[m] = exp(i*max(d[d<0])/sum_smaller)
    i=i+1
  }
  return(z) 
}

integration <- function(x,yf){
  n <- length(x)
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  return(integral)
}

# PAFFT
PAFFT <- function(spectra,reference,segSize,shift){
  lags <- alignedSpectrum <- matrix(0,nrow(spectra),ncol(spectra))
  for (i in 1:nrow(spectra)){
    startpos <- 1
    aligned <- c()
    shifts <- rep(0,ncol(spectra))
    while (startpos <= ncol(spectra)){
      endpos <- startpos+(segSize*2)
      if (endpos>=ncol(spectra)){
        samseg <- spectra[i,startpos:ncol(spectra)]
        refseg <- reference[startpos:ncol(spectra)]
      }else{
        samseg <- spectra[i,(startpos+segSize):(endpos-1)]
        refseg <- reference[(startpos+segSize):(endpos-1)]
        minpos <- findMin(samseg,refseg)
        endpos <- as.numeric(startpos+minpos+segSize)
        samseg <- spectra[i,startpos:endpos]
        refseg <- reference[startpos:endpos]
      }
      lag <- FFTcorr(samseg,refseg,shift)
      samseg <- as.numeric(samseg)
      shifts[startpos:(startpos+length(refseg)-1)] <- lag
      aligned <- c(aligned,as.numeric(move(samseg,lag)))
      startpos <- endpos+1
    }
    alignedSpectrum[i,] <- aligned
    lags[i,] <- shifts
  }
  return(list(alignedSpectrum=alignedSpectrum,lags=lags))
}

FFTcorr <- function(spectrum, target, shift){
  spectrum <- t(as.matrix(spectrum))
  target <- t(as.matrix(target))
  M <- ncol(target)
  diff <- 1000000
  for (i in 1:20){
    curdiff <- ((2^i)-M)
    if (curdiff>0&curdiff<diff){diff <- curdiff}
  }
  target <- cbind(target,t(rep(0,diff)))
  spectrum <- cbind(spectrum,t(rep(0,diff)))
  M <- M+diff
  X <- fft(as.numeric(target))
  Y <- fft(as.numeric(spectrum))
  R <- X*Conj(Y)
  R <- R/M
  rev <- fft(R,inverse=T)/length(rev)
  vals <- Re(rev)
  maxpos <- 1
  maxi <- -1
  if (M<shift){shift <- M}
  for (i in 1:shift){
    if (vals[i] > maxi){
      maxi = vals[i]
      maxpos = i
    }
    if (vals[length(vals)-i+1] > maxi){
      maxi = vals[length(vals)-i+1];
      maxpos = length(vals)-i+1;
    }
  }
  if (maxi < 0.1){lag <- 0
  return(lag)}
  if (maxpos > length(vals)/2){
    lag = maxpos-length(vals)-1
    return(lag)
  }else{lag <- maxpos-1
  return(lag)}
}

move <- function(seg, lag){
  if (lag == 0 || lag >= length(seg)){
    movedSeg <- seg
    return(movedSeg)
  }
  if (lag > 0){
    ins <- rep(1,lag)*seg[1]
    movedSeg <- c(ins,seg[1:(length(seg)-lag)])
    return(movedSeg)
  }
  if (lag < 0){
    lag <- abs(lag);
    ins <- rep(1,lag)*seg[length(seg)]
    movedSeg <- c(seg[(lag+1):length(seg)],ins)
    return(movedSeg)
  }
}

findMin <- function(samseg,refseg){
  Is <- order(samseg)
  Ir <- order(refseg)
  minposA <- c()
  minInt <- c()
  minpos <- NA
  for (i in 1:round(length(Is)/20)){
    for (j in 1:round(length(Is)/20)){
      if (Ir[j]==Is[i]){minpos <- Is[i]
      return(minpos)}
    }
  }
  if (is.na(minpos)){minpos <- Is[1]}
  return(minpos)
}
# End PAFFT