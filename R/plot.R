plotEICs <- function(eics, peaks=NULL, rtrange=NULL) {
  p <- plot_ly() %>%
    layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
           yaxis = list(title = 'Intensity'))
  
  for (f in 1:length(eics)) {
    p <- add_trace(p, x = eics[[f]]$rt, y = eics[[f]]$intensity, mode='line', name = names(eics)[f] )
  }
  
  if (!is.null(peaks)){
    eic_mat <- do.call(rbind, lapply(eics, function(e) e$intensity))
    eic <- list(rt=eics[[1]]$rt, intensity=apply(eic_mat, 2, max))
    p <- add_markers(p, x = eic$rt[peaks$Index$Position], y = eic$intensity[peaks$Index$Position], name = 'peak position', color = I('red'), marker = list(size = 5))
    p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = eic$intensity[c(peaks$Index$Start, peaks$Index$End)], name = 'peak bound', color = I('blue'), marker = list(size = 5))
  }
  
  if (!is.null(rtrange)){
    eic_mat <- do.call(rbind, lapply(eics, function(e) e$intensity))
    eic <- list(rt=eics[[1]]$rt, intensity=apply(eic_mat, 2, max))
    xx <- findInterval(rtrange,eic$rt)
    xx[xx<1] <- 1
    yy <- eic$intensity[xx]
    p <- add_markers(p, x = rtrange, y = yy, name = 'peak bound', color = I('blue'), marker = list(size = 5))
  }
  
  p
}

plotStackBar <- function(SampleAreas) {
  SampleNames <- colnames(SampleAreas)
  SampleRelAreas <- t(t(SampleAreas)/colSums(SampleAreas))
    
  p <- plot_ly(x = SampleNames, y = SampleRelAreas[1,], name=rownames(SampleRelAreas)[1], type = 'bar') %>%
    layout(barmode = 'stack',
           xaxis = list(title = 'Sample Name'),
           yaxis = list(title = 'Relative Abundance'))
  
  for (i in 2:nrow(SampleAreas)) {
    p <- add_trace(p, x = SampleNames, y = SampleRelAreas[i,], name=rownames(SampleRelAreas)[i])
  }
  
  p
}