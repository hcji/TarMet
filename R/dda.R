getMS2.DDA <- function(rawDIADataset, target.mz, target.rt, dda.mztol=0.01, dda.rttol=5, abund.Th=0.01) {
  MS2 <- lapply(rawDIADataset, function(this){
    pmz <- {
      res <- as.numeric(unlist(strsplit(names(this), 'precursorMZ:')))
      res[!is.na(res)]
    }
    prt <- vapply(this, function(s) s$times[1], 0)
    
    wh <- which(abs(pmz - target.mz) < dda.mztol & abs(prt - target.rt) < dda.rttol)
    if (length(wh) < 1) {
      return(NULL)
    } else {
      wh <- wh[which.min(abs(pmz[wh] - target.mz))]
    }
    
    ms2 <- this[[wh]]$peaks[[1]]
    ms2 <- ms2[ms2[,2]>abund.Th*max(ms2[,2]), ]
    ms2[,2] <- ms2[,2]/max(ms2[,2])
    return(ms2)
  })
  
  names(MS2) <- names(rawDIADataset)
  return(MS2)
}