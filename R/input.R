getSampleName <- function(files){
  Names <- sapply(files, function(f){
    Name <- strsplit(f,'/')[[1]]
    Name[length(Name)]
  })
  return(as.character(Names))
}

getCandidateFormula <- function(mz, window, adduct) {
  wh <- which(adduct==adducts$Name)
  charge <- abs(adducts$Charge[wh])
  mass <- mz * charge - adducts$Mass[wh]
  formulas <- generate.formula(mass, window=window, validation=TRUE, charge=charge)
  out <- sapply(formulas, function(f) f@string)
  return(out)
}

getTargets <- function(config, define) {
  info <- read.csv(config)
  if (define=='mass-to-charge'){
    return(info$mz)
  } else {
    return(as.character(info$formula))
  }
}

getElements <- function(formula) {
  formula <- rcdk::get.formula(formula)
  return(formula@isotopes[,1])
}

getElementNum <- function(formula, element){
  formula <- rcdk::get.formula(formula)
  return(as.numeric(formula@isotopes[formula@isotopes[,1] == element, 2]))
}