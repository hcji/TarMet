runTarMet <- function(){
  appdir <- system.file('app', package = 'TarMet')
  runApp(appdir, display.mode = 'normal')
}

LoadData <- function(filename)
{
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

getIsoPat <- function(formula, d, threshold, resolution){
  data("isotopes", package = "enviPat")
  data("adducts", package = "enviPat")
  if (!is.numeric(d)){d <- which(d == adducts$Name)}
  checked <- check_chemform(isotopes, formula)
  if (d > 0){
    if (as.numeric(adducts[d, 4]) == 2){
      chemform <- mergeform(checked[,2], checked[,2])
    } else if (as.numeric(adducts[d, 4]) == 3){
      chemform <- mergeform(checked[,2], checked[,2], checked[,2])
    } else {
      chemform <- checked[, 2]
    }
    if (adducts[d, 7] != 'FALSE'){
      chemform <- mergeform(chemform, adducts[d, 7])
    } 
    if (adducts[d, 8] != 'FALSE'){
      chemform <- subform(chemform, adducts[d, 8])
    }
  }
  checked <- check_chemform(isotopes, chemform)
  pattern <- isowrap(isotopes, checked, resmass=FALSE, threshold = threshold, resolution = resolution)[[1]]
  return(as.data.frame(pattern))
}

getIsoEIC.mz <- function(raw, fmz, nmax = 4, ppm = 50, rtrange = c(0, Inf), charge = 1,
                         C13 = 3, H2 = 0, O18 = 0, N15 = 0, S34 = 0){
  
  mass <- (0:C13 * 1.003355)
  mass <- unlist(lapply(0:H2 * 1.006277, function(s) s+mass))
  mass <- unlist(lapply(0:O18 * 2.004245, function(s) s+mass))
  mass <- unlist(lapply(0:N15 * 0.997035, function(s) s+mass))
  mass <- unlist(lapply(0:S34 * 1.995796, function(s) s+mass))
  mass <- mass[mass < (nmax + 0.1)]
  mass <- sort(mass)

  mzs <- (fmz + mass)/charge
  remain <- diff(mzs)/mzs[-1]*10^6 > 0.5 * ppm
  mzs <- mzs[remain]
  
  mzranges <- do.call(rbind, lapply(mzs, function(mz){
    c(mz * (1 - ppm/2/10^6), mz * (1 + ppm/2/10^6))
  }))
  eics <- apply(mzranges, 1, function(mzrange){
    getEIC(raw, rtrange, mzrange)
  })
  
  return(list(mzs = mzranges, eics = eics))
}

getIsoEIC.formula <- function(raw, formula, adduct = 'M+H', ppm = 50, rtrange = c(0, Inf), threshold = 0.01, resolution = 50000){

  data("adducts", package = "enviPat")
  adduct <- which(adduct == adducts$Name)
  
  pattern <- getIsoPat(formula, adduct, threshold, resolution)
  mzs <- pattern[,1]
  ppm <- min(ppm, 0.5*(diff(pattern[,1])/pattern[-1,1]*10^6))

  mzranges <- do.call(rbind, lapply(mzs, function(mz){
    c(mz * (1 - ppm/2/10^6), mz * (1 + ppm/2/10^6)) / abs(adducts$Charge[adduct])
  }))
  
  eics <- apply(mzranges, 1, function(mzrange){
    getEIC(raw, rtrange, mzrange)
  })
  
  return(list(mzs = mzranges, pattern = pattern, eics = eics))
}

getIsoPeaks <- function(eics, SNR.Th = 4, peakScaleRange = 5, peakThr = 0, userDefTime = NULL){
  peakScaleRange <- round(peakScaleRange/ mean(diff(eics$eics[[1]]$rt)))
  ind <- which.max(sapply(eics$eics, function(s){
    sum(s$intensity)
  }))
  eic <- eics$eics[[ind]]
  if (sum(eic$intensity)==0) return(NULL)
  MajorPeaks <- peakDetectionCWT(eic$intensity, scales = c(1, seq(2, round(0.05*length(eic$intensity)), 2)), SNR.Th = SNR.Th, peakScaleRange = peakScaleRange, peakThr = peakThr)
  if (length((MajorPeaks$majorPeakInfo$peakIndex)) < 1) {return(NULL)}
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

getArea <- function(eics, rtmin, rtmax, intensity = FALSE){
  areas <- sapply(eics$eics, function(eic){
    Start <- findInterval(rtmin, eic$rt)+1
    End <- findInterval(rtmax, eic$rt)
    if (intensity == TRUE){
      Area <- max(eic$intensity[Start:End])
    } else {
      Area <- integration(eic$rt[Start:End], eic$intensity[Start:End])
    }
    
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

plotEICs <- function(eics, peaks = NULL, Names = NULL, rt = NULL) {
  if (!is.null(rt)){
    eics <- lapply(eics, function(eic){
      ids <- eic$rt>rt[1] & eic$rt<rt[2]
      eic$rt <- eic$rt[ids]
      eic$intensity <- eic$intensity[ids]
      eic
    })
  }
  p <- plot_ly() %>%
    layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
           yaxis = list(title = 'Intensity'))
  if (is.null(Names)){
    Names = paste('Mz: ',round(eics$mzs[,1],4), ' - ', round(eics$mzs[,2],4))
  }
  for (f in 1: length(eics$eics)) {
    p <- add_trace(p, x = eics$eics[[f]]$rt, y = eics$eics[[f]]$intensity, mode='line', name = paste(Names[f]))
  }
  if (!is.null(peaks)){
    eic <- eics$eics[[1]]
    p <- add_markers(p, x = eic$rt[peaks$Index$Position], y = eic$intensity[peaks$Index$Position], name = 'peak position', color = I('red'), marker = list(size = 5))
    p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = eic$intensity[c(peaks$Index$Start, peaks$Index$End)], name = 'peak bound', color = I('blue'), marker = list(size = 5))
  }
  return(p)
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
# End PAFFT

# clone from baselineWavelet

cwt <- function(ms, scales=1, wavelet='mexh') {
  if (wavelet == 'mexh') {
    psi_xval <- seq(-8, 8, length=1024)
    psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *exp(-psi_xval^2/2)
    #plot(psi_xval, psi)
  } else if (wavelet=='haar') {
    psi_xval <- seq(0,1,length=1024)
    psi <- c(0,rep(1,511),rep(-1,511),0)
  }else if (is.matrix(wavelet)) {
    if (nrow(wavelet) == 2) {
      psi_xval <- wavelet[1,]
      psi <- wavelet[2,]
    } else if (ncol(wavelet) == 2) {
      psi_xval <- wavelet[,1]
      psi <- wavelet[,2]
    } else {
      stop('Unsupported wavelet format!')
    }
  } else {
    stop('Unsupported wavelet!')
  }
  
  oldLen <- length(ms)
  ms <- extendNBase(ms, nLevel=NULL, base=2)
  len <- length(ms)
  nbscales <- length(scales)
  wCoefs <- NULL
  
  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax  <- psi_xval[length(psi_xval)]
  for (i in 1:length(scales)) {
    scale.i <- scales[i]
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
    if (length(j) == 1)		j <- c(1, 1)
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
    wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
    ## Shift the position with half wavelet width
    wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
    wCoefs <- cbind(wCoefs, wCoefs.i)
  }
  if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
  colnames(wCoefs) <- scales
  wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
  return(wCoefs)
}

WhittakerSmooth <- function(x,w,lambda,differences=1) {
  x=as.vector(x)
  L=length(x)
  E=spMatrix(L,L,i=seq(1,L),j=seq(1,L),rep(1,L))
  D=as(diff(E,1,differences),"dgCMatrix")
  W=as(spMatrix(L,L,i=seq(1,L),j=seq(1,L),w),"dgCMatrix")
  background=solve((W+lambda*t(D)%*%D),w*x);
  return(as.vector(background))
}

widthEstimationCWT <- function(x,majorPeakInfo) {
  
  wCoefs_haar <- cwt(x, 1:max(majorPeakInfo$peakScale), wavelet='haar')
  peakIndex <- majorPeakInfo$peakIndex
  peakScale <- majorPeakInfo$peakScale[findInterval(majorPeakInfo$peakIndex,majorPeakInfo$allPeakIndex)]
  
  peakWidth <- list()
  peakWidth[["peakIndex"]]= majorPeakInfo$peakIndex
  
  for(i in 1:length(peakIndex)){
    peakIndex.i=peakIndex[i]
    peakScale.i=peakScale[i]
    wCoefs_haar.i=wCoefs_haar[,peakScale.i]
    wCoefs_haar.i.abs=abs(wCoefs_haar.i)
    localmax=localMaximum(-wCoefs_haar.i.abs,winSize=5)      
    #    localmax=localmax & abs(wCoefs_haar.i)<(mean(wCoefs_haar.i.abs[localmax==1])+0.5*sd(wCoefs_haar.i.abs[localmax==1]))
    localmax=as.numeric(localmax)
    localmax[peakIndex]=0
    localmax[(peakIndex.i-peakScale.i/2+1):(peakIndex.i+peakScale.i/2-1)]=0
    
    Lef =0
    Rig =0
    
    peakScale.i.3=2*peakScale.i
    
    
    if(i==1){
      maxIndexL=max(c((peakIndex.i-peakScale.i.3),1))
    }else{
      maxIndexL=max(c((peakIndex.i-peakScale.i.3),peakIndex[i-1]))
    }
    
    if(i==length(peakIndex)){
      minIndexR=min(c((peakIndex.i+peakScale.i.3),length(localmax)))
    } else{
      minIndexR=min(c((peakIndex.i+peakScale.i.3),peakIndex[i+1]))
    }
    ignoreL=1:maxIndexL
    ignoreR=minIndexR:length(localmax)        
    localmax[c(ignoreL,ignoreR)]=0
    localmax[c(peakIndex.i,(peakIndex.i-(peakScale.i/2)):(peakIndex.i+(peakScale.i/2)))]=0
    bi=which(localmax==1)
    
    biLeft=bi[bi<peakIndex.i]
    useL= maxIndexL:peakIndex.i
    minIndexLeft=useL[which(min(x[useL])==x[useL])]
    
    if(length(biLeft)==0){
      Lef=minIndexLeft
    }else{
      minbaselineIndexLeft=biLeft[which(min(x[biLeft])==x[biLeft])]
      if(minIndexLeft>=(peakIndex.i-peakScale.i/2+1)){
        Lef=minbaselineIndexLeft
      }else{
        Lef=max(c(minIndexLeft,minbaselineIndexLeft))
      }   
    }
    
    biRight=bi[bi>peakIndex.i]
    useR= peakIndex.i:minIndexR
    minIndexRight=useR[which(min(x[useR])==x[useR])]
    
    if(length(biRight)==0){
      Rig= minIndexRight
    }else{
      minbaselineIndexRight=biRight[which(min(x[biRight])==x[biRight])] 
      if(minIndexRight<=(peakIndex.i+peakScale.i/2-1)){
        Rig=minbaselineIndexRight
      }else{
        Rig=min(c(minIndexRight,minbaselineIndexRight))
      }      
    }
    
    peakWidth[[paste(peakIndex.i)]]=Lef:Rig
    
  }
  
  return(peakWidth)		
}

extendNBase <-  function(x, nLevel=1, base=2, ...) {
  if (!is.matrix(x)) x <- matrix(x, ncol=1)	
  
  nR <- nrow(x)
  if (is.null(nLevel)) {
    nR1 <- nextn(nR, base)		
  } else {
    nR1 <- ceiling(nR / base^nLevel) * base^nLevel		
  }
  if (nR != nR1) {
    x <- extendLength(x, addLength=nR1-nR, ...)
  }
  
  return(x)
}

extendLength <- function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both')) {
  if (is.null(addLength)) stop('Please provide the length to be added!')
  if (!is.matrix(x)) x <- matrix(x, ncol=1)	
  method <- match.arg(method)
  direction <- match.arg(direction)
  
  nR <- nrow(x)
  nR1 <- nR + addLength
  if (direction == 'both') {
    left <- right <- addLength
  } else if (direction == 'right') {
    left <- 0
    right <- addLength
  } else if (direction == 'left') {
    left <- addLength
    right <- 0
  }
  
  if (right > 0) {
    x <- switch(method,
                reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
                open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
                circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
  }
  
  if (left > 0) {
    x <- switch(method,
                reflection =rbind(x[addLength:1, , drop=FALSE], x),
                open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
                circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
  }
  if (ncol(x) == 1)  x <- as.vector(x)
  
  return(x)
}
