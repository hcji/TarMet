getQuantifiedResult <- function(SamplesEICs, targetRtLeft, targetRtRight){
  res <- lapply(SamplesEICs, function(EICs){
    round(getArea(EICs, targetRtLeft, targetRtRight), 2)
  })
  res <- do.call(cbind, res)
  return(res)
}