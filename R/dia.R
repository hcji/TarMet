loadSWATH <- function(datapath) {
  # check format
  splitname <- strsplit(datapath,"\\.")[[1]]
  if (tolower(splitname[length(splitname)]) != "mzml"){
    msobj <- openMSfile(datapath, backend="pwiz")
  } else if (tolower(splitname[length(splitname)]) != "mzxml"){
    msobj <- openMSfile(datapath, backend="Ramp")
  } else {
    stop('Only .mzML and .mzXML format is supported')
  }
  
  # open ms file
  msobj <- openMSfile(datapath, backend="pwiz")
  peakInfo <- peaks(msobj)
  headerInfo <- header(msobj)
  
  # ms2 data
  precursorMZ <- unique(headerInfo$precursorMZ)
  precursorMZ <- precursorMZ[precursorMZ!=0]
  dataMS2 <- lapply(precursorMZ, function(pmz){
    this <- headerInfo$precursorMZ == pmz
    peakInfo.this <- peakInfo[this]
    
    peakInfo.this <- lapply(peakInfo.this, function(spectrum) {
      keep <- spectrum[,2] > 1e-6
      output <- as.data.frame(spectrum[keep,,drop = FALSE])
      colnames(output) <- c('mz','intensity')
      return(output)
    })
    
    data.this <- list(
      peaks = peakInfo.this,
      times = headerInfo$retentionTime[this]
    )
  })
  names(dataMS2) <- paste('precursorMZ:', precursorMZ, sep='')
  
  return(dataMS2)
}

getEIC.SWATH <- function(rawDIADataset, targetMz, targetPeaks, resolution, ppm) {
  peaks <- targetPeaks$PeakInfo
  diaEICs <- list()
  for (i in 1:nrow(peaks)){
    this.start <- as.numeric(peaks$Start[i])
    this.end <- as.numeric(peaks$End[i])
    this.position <- as.numeric(peaks$Position[i])
    
    this.data <- lapply(rawDIADataset, function(this.SWATH){
      pmz <- {
        res <- as.numeric(unlist(strsplit(names(this.SWATH), 'precursorMZ:')))
        res[!is.na(res)]
      }
      wh <- which.min(abs(pmz-targetMz))
      this.SWATH[[wh]]
    })
    
    this.eics <- lapply(this.data, function(this){
      scan <- which.min(abs(this$times-this.position))
      mz <- this$peaks[[scan]]$mz
      mzranges <- getMzRanges(mz, resolution, ppm)
      rtranges <- c(this.start, this.end)
      getMzEICs(this, mzranges, rtranges, baseline=FALSE, smooth=FALSE)
    })
    
    diaEICs <- c(diaEICs, list(this.eics))
  }
  
  names(diaEICs) <- peaks$Name
  return(diaEICs)
}

getMS2.SWATH <- function(targetEICs, targetPeaks, diaEICs, msCorr.Th=0.8){
  peaks <- targetPeaks$PeakInfo
  ms2 <- list()
  for (i in 1:nrow(peaks)){
    this.start <- as.numeric(peaks$Start[i])
    this.end <- as.numeric(peaks$End[i])
    this.position <- as.numeric(peaks$Position[i])
    
    this.eic.ms1 <- lapply(targetEICs, function(eics){
      wh <- eics[[1]]$rt >= this.start & eics[[1]]$rt <= this.end
      list(rt=eics[[1]]$rt[wh], intensity=eics[[1]]$intensity[wh])
    })
    this.eics.ms2 <- diaEICs[[i]]
    
    this.ms2 <- lapply(seq_along(this.eic.ms1), function(s){
      eic1 <- this.eic.ms1[[s]]
      eic2 <- this.eics.ms2[[s]]
      
      scores <- scoreEICs(eic1, eic2)
      mass <- sapply(names(eic2), function(mzr){
        mean(as.numeric(strsplit(mzr, '-')[[1]]))
      })
      intensity <- sapply(eic2, function(eic){
        wh <- which.min(abs(eic$rt - this.position))
        eic$intensity[wh]
      })
      
      keep <- which(scores>=msCorr.Th)
      if (length(keep)>1){
        data.frame(mz=mass[keep], intensity=intensity[keep], scores=scores[keep])
      } else {
        NULL
      }
    })
    ms2 <- c(ms2, this.ms2)
  }
  return(ms2)
}

scoreEICs <- function(eic1, eic2, scales=1:24, points=500){
  freq <- (max(eic1$rt)-min(eic1$rt)) / points
  x <- seq(min(eic1$rt), max(eic1$rt), freq)
  scores <- sapply(eic2, function(this){
    y1 <- approx(x=eic1$rt, y=eic1$intensity, xout=x, rule=2)$y
    y2 <- approx(x=this$rt, y=this$intensity, xout=x, rule=2)$y
    ym1 <- conv(y1, scales)
    ym2 <- conv(y2, scales)
    matcor(ym1, ym2)
  })
  return(scores)
}

getScores.SWATH <- function(formula, MS2, tarID, msDB, ppm=10, adduct='M+H', typeDB='experimental', eval='median'){
  scores <- lapply(MS2, function(ms2){
    getMatchScore(formula, ms2, tarID, msDB, ppm, adduct, typeDB)
  })
  type <- scores[[1]]$type
  matching <- sapply(scores, function(s){s$scores['matching']})
  correlation <- sapply(scores, function(s){s$scores['correlation']})
  if (eval=='mean'){
    matching <- mean(matching)
    correlation <- mean(correlation)
  } else {
    matching <- median(matching)
    correlation <- median(correlation)
  }
  return(list(type=type, matching=matching, correlation=correlation))
}
