getIsotopicPeaks <- function(eics, SNR.Th = 4, peakScaleRange = 5, peakThr = 0, theoretical = NULL){
  peakScaleRange <- round(peakScaleRange/ mean(diff(eics[[1]]$rt)))
  eic_mat <- do.call(rbind, lapply(eics, function(e) e$intensity))
  eic <- list(rt=eics[[1]]$rt, intensity=apply(eic_mat, 2, max))
  if (sum(eic$intensity)==0) return(NULL)
  MajorPeaks <- try(
    peakDetectionCWT(eic$intensity, scales = c(1, seq(2, round(0.05*length(eic$intensity)), 2)), SNR.Th = SNR.Th, peakScaleRange = peakScaleRange, peakThr = peakThr)
  )
  if (class(MajorPeaks)=='try-error'){
    return(NULL)
  }
  if (length((MajorPeaks$majorPeakInfo$peakIndex)) < 1) {return(NULL)}
  PeakWidths <- widthEstimationCWT(eic$intensity, MajorPeaks$majorPeakInfo)
  
  Position <- unique(MajorPeaks$majorPeakInfo$peakIndex)
  Start <- sapply(2:length(PeakWidths), function(s){
    PeakWidths[[s]][1]
  })
  End <- sapply(2:length(PeakWidths), function(s){
    PeakWidths[[s]][length(PeakWidths[[s]])]
  })
  
  Area <- do.call(rbind, lapply(eics, function(eic){
    sapply(1:length(Position), function(p){
      integration(eic$rt[Start[p]:End[p]], eic$intensity[Start[p]:End[p]] )
    })
  }))
  if (!is.null(theoretical)){
    Similarity <- sapply(1:ncol(Area), function(s){
      round(getIsotopicSimilarity(Area[,s], theoretical),3)
    })
  } else {
    Similarity <- NA
  }
  PeakArea <- as.data.frame(round(Area,2))
  colnames(PeakArea) <- c(paste('Peak', 1:length(Position)))
  PeakInfo <- data.frame(Name = paste('Peak', 1:length(Position)), 
                         Position = eic$rt[Position], 
                         Start = eic$rt[Start], 
                         End = eic$rt[End], 
                         Similarity = Similarity)
  Index <- data.frame(Position=Position, Start=Start, End=End)
  
  return(list(PeakInfo = PeakInfo, PeakArea = PeakArea, Index = Index))
}

getArea <- function(eics, rtmin, rtmax){
  areas <- sapply(eics, function(eic){
    Start <- findInterval(rtmin, eic$rt)+1
    End <- findInterval(rtmax, eic$rt)
    Area <- integration(eic$rt[Start:End], eic$intensity[Start:End])
    Area
  })
  
  return(areas)
}

integration <- function(x,yf){
  n <- length(x)
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  return(integral)
}

getIsotopicSimilarity <- function(iso_int1, iso_int2) {
  iso_int1 <- (iso_int1/sum(iso_int1)) [iso_int2!=0]
  iso_int2 <- (iso_int2/sum(iso_int2)) [iso_int2!=0]
  score_i <- sapply(seq_along(iso_int1), function(s){
    (1-abs(iso_int1[s]-iso_int2[s])/iso_int2[s])*iso_int2[s]
  })
  score <- sum(score_i)
  
  return(score)
}

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