getAlignedEICs <- function(eics, ind, shift, segment) {
  if (length(eics)==1 | shift<0 | segment<1){
    return(eics)
  }
  
  rtout <- eics[[ind]][[1]]$rt
  references <- lapply(eics[[ind]], function(eic) eic$intensity)
  
  segSize <- round(segment/mean(diff(rtout)))
  Shift <- round(shift/mean(diff(rtout)))
  
  for (i in seq_along(eics)[-ind]) {
    for (j in seq_along(references)){
      spectra <- eics[[i]][[j]]
      spectra <- approx(x=spectra$rt, y=spectra$intensity, xout=rtout, rule = 2)$y
      aligned <- PAFFT(t(spectra), t(references[[j]]), segSize, Shift)$alignedSpectrum
      eics[[i]][[j]]$rt <- rtout
      eics[[i]][[j]]$intensity <- as.numeric(aligned)
    }
  }
  
  return(eics)
}

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
  rev <- fft(R,inverse=TRUE)/length(rev)
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