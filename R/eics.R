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

getEIC <- function(raw, rtrange=c(0,Inf), mzrange){
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

getIsotopicPattern <- function(formula, adduct, threshold, resolution){
  data("isotopes", package = "enviPat")
  data("adducts", package = "enviPat")
  d <- adduct
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

getMzWithTracer <- function(formula, adduct, element = c('C'), isotope = c('13C'), number = c(3)){
  data("isotopes", package = "enviPat")
  data("adducts", package = "enviPat")
  d <- adduct
  if (!is.numeric(d)){d <- which(d == adducts$Name)}
  
  mass <- rcdk::get.formula(formula)@mass
  mass <- (mass + adducts$Mass[d])
  
  addmass <- 0
  if (!is.null(element) && length(element)==length(number) && length(element) == length(isotope)) {
    for (i in seq_along(element)) {
      m1 <- isotopes$mass[element[i] == isotopes$element][1]
      m2 <- isotopes$mass[isotope[i] == isotopes$isotope][1]
      addmass <- unlist(lapply(0:number[i] * (m2 - m1), function(s) s+addmass))
    }
  }
  charge <- abs(adducts$Charge[d])
  mzs <- (mass+addmass)/charge
  
  return(mzs)
}

getMzRanges <- function(mzs, resolution){
  mzranges <- do.call(rbind, lapply(mzs, function(mz){
    delta <- mz/resolution
    data.frame(mzmin=mz-delta, mzmax=mz+delta)
  }))
  return(mzranges)
}

getMzEICs <- function(raw, mzranges, baseline=TRUE, smooth=FALSE){
  out <- lapply(1:nrow(mzranges), function(s){
    eic <- getEIC(raw, mzrange=as.numeric(mzranges[s,]))
    if (smooth) {
      eic$intensity <- smooth.spline(eic$intensity)$y
    }
    if (baseline){
      eic$intensity <- eic$intensity - airPLS(eic$intensity)
    }
    eic
  })
  names(out) <- paste(round(mzranges[,1],4), round(mzranges[,2],4), sep = '-')
  return(out)
}
