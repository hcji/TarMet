getExperMS <- function(tarID, msDB, adduct='M+H', typeDB='experimental'){
  tarID <- as.character(tarID)
  if (typeDB=='experimental'){
    wh <- which(as.character(msDB$id)==tarID & msDB$adduct==adduct)
    if (length(wh>2)){
      return(msDB[wh,c('ProductMz', 'LibraryIntensity')])
    } else {
      return(NULL)
    }
  }
}

getMatchScore <- function(ms2, tarID, msDB, ppm=50, adduct='M+H', typeDB='experimental') {
  if (typeDB == 'experimental') {
    ms2_std <- getExperMS(tarID, msDB, adduct, typeDB)
    ms2_match <- lapply(1:nrow(ms2_std), function(i){
      mz.std <- ms2_std[i,1]
      mz.diff <- abs(ms2[,'mz']-mz.std)/mz.std*10^6
      if (min(mz.diff)>ppm){
        mz <- mz.std
        intensity <- 0
        corr <- NA
      } else {
        wh <- which.min(mz.diff)
        mz <- ms2[wh,'mz']
        intensity <- ms2[wh,'intensity']
        corr <- ms2[wh,'scores']
      }
      c(mz=mz, intensity=intensity, corr=corr)
    })
    ms2_match <- do.call(rbind, ms2_match)
    ms2_match[,2] <- ms2_match[,2]/max(ms2_match[,2])
    matching.int <- (1-abs(ms2_std[,2]-ms2_match[,2]) / apply(rbind(ms2_match[,2], ms2_std[,2]),2,max))
    matching.mz <- (1-abs(ms2_std[,1]-ms2_match[,1]))
    matching <- sum(matching.mz^1.84 * matching.int^0.59* ms2_std[,2]) / sum(ms2_std[,2])
    correlation <- mean(ms2_match[,3], na.rm = TRUE)
    if (is.nan(correlation)) {correlation <- 0}
    type <- 'experimental'
    score <- c(matching=matching, correlation=correlation)
  }
  return(list(type=type, stdMS = ms2_std, score=score))
}