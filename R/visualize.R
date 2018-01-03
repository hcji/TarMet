plotEICs <- function(eics, peaks = NULL, Names = NULL) {
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