plotEICs <- function(eics, peaks) {
  p <- plot_ly() %>%
    layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
           yaxis = list(title = 'Intensity'))
  
  for (f in 1: length(eics)) {
    p <- add_trace(p, x = eics[[f]]$rt, y = eics[[f]]$intensity, mode='line', name = names(eics)[f] )
  }
  
  if (!is.null(peaks)){
    eic_mat <- do.call(rbind, lapply(eics, function(e) e$intensity))
    eic <- list(rt=eics[[1]]$rt, intensity=apply(eic_mat, 2, max))
    p <- add_markers(p, x = eic$rt[peaks$Index$Position], y = eic$intensity[peaks$Index$Position], name = 'peak position', color = I('red'), marker = list(size = 5))
    p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = eic$intensity[c(peaks$Index$Start, peaks$Index$End)], name = 'peak bound', color = I('blue'), marker = list(size = 5))
  }
  
  p
}