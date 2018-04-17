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

getMatchScore <- function(formula, ms2, tarID, msDB, adduct='M+H', typeDB='experimental') {
  if (typeDB == 'experimental') {
    experMS <- getExperMS(tarID, msDB, adduct='M+H', typeDB='experimental')
  } else {
    typeDB <- 'in-silicon'
  }
  
  if (typeDB=='in-silicon'){
    settingsObject<-list()
    settingsObject[["DatabaseSearchRelativeMassDeviation"]] <- 10
    settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]] <- 0.001
    settingsObject[["FragmentPeakMatchRelativeMassDeviation"]] <- 5.0
    settingsObject[["MetFragDatabaseType"]] <- 'PubChem'
    settingsObject[["PeakList"]] <- ms2
    settingsObject[["NeutralPrecursorMolecularFormula"]] <- formula
    settingsObject[["PrecursorCompoundIDs"]] <- as.character(tarID)
    
    scored.candidates<-run.metfrag(settingsObject)
    score <- scored.candidates$Score[1]
  } else {
    
  }
  
}