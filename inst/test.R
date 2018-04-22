input <- list()
input$type <- 'data independent analysis'
input$files$datapath <- list.files('E:/dataset/MetDIA_demo/data', full.names = TRUE)
input$files$name <- list.files('E:/dataset/MetDIA_demo/data', full.names = TRUE)

# define samples
sampleNames <- {
  req(input$files)
  getSampleName(input$files$name)
}

rawDataset <- {
  res <- list()
  for (i in seq_along(sampleNames)){
    res[[i]] <- LoadData(input$files$datapath[i])
  }
  names(res) <- sampleNames
  res
}

rawDIADataset <- {
  if (input$type == 'data independent analysis') {
    res <- list()
    for (i in seq_along(sampleNames)){
      res[[i]] <- loadSWATH(input$files$datapath[i])
    }
    names(res) <- sampleNames
    res
  }
}

input$formula <- 'C7H15NO3'
input$adduct <- 'M+H'
input$ppm <- 20
input$rtmin <- 0
input$rtmax <- Inf
input$snr.th <- 5
input$scale.th <- 5
input$int.th <- 0
input$resolution <- 50000
input$threshold <- 0.01
sampleInd <- 1

formula <- gsub(" ", "", as.character(input$formula))
pattern <- getIsotopicPattern(formula, input$adduct, input$threshold, input$resolution)

targetMzRanges <- {
  if (input$type=='isotopic tracer'){
    mzs <- getMzWithTracer(formula, input$adduct, input$tracer_element, input$tracer_isotope, input$tracer_number)
  } else {
    mzs <- pattern[,1]
  }
  getMzRanges(mzs, resolution=input$resolution, ppm=input$ppm)
}

targetEICs <-{
  res <- list()
  rtranges <- c(input$rtmin, input$rtmax)
  for (i in seq_along(sampleNames)){
    res[[i]] <- getMzEICs(rawDataset[[i]], rtranges=rtranges, mzranges=targetMzRanges)
  }
  names(res) <- sampleNames
  res
}

targetPeaks <- {
  if (input$type!='isotopic tracer'){
    theoretical <- as.numeric(pattern[,2])
  } else {
    theoretical <- NULL
  }
  getIsotopicPeaks(targetEICs[[sampleInd]], SNR.Th=input$snr.th, peakScaleRange=input$scale.th, peakThr=input$int.th, theoretical=theoretical)
}

whichPeak <- {
  if (input$type!='isotopic tracer'){
    which.max(targetPeaks$PeakInfo$Similarity)
  } else {
    which.max(colSums(targetPeaks$PeakArea))
  }
}

input$targetRtPosition <- targetPeaks$PeakInfo$Position[whichPeak]
input$targetRtLeft <- targetPeaks$PeakInfo$Start[whichPeak]
input$targetRtRight <- targetPeaks$PeakInfo$End[whichPeak]
input$msCorr.Th <- 0.8
input$tarID <- "HMDB00062"
input$msDB$datapath <- 'E:/project/TarMet/inst/data/DemoDatabase.csv'

userInfo <- {
  Position <- input$targetRtPosition
  Start <- input$targetRtLeft
  End <- input$targetRtRight
  data.frame(Name='User', Position=Position, Start=Start, End=End)
}

outputPeakInfo <- {
  rbind(targetPeaks$PeakInfo[,1:ncol(userInfo)], userInfo)
}

msDB <- read.csv(input$msDB$datapath)

diaEICs <- {
  targetMz <- mean(as.numeric(targetMzRanges[1,]))
  getEIC.SWATH(rawDIADataset, targetMz, outputPeakInfo, input$resolution, input$ppm)
}

MS2 <- getMS2.SWATH(targetEICs, outputPeakInfo, diaEICs, input$msCorr.Th)
diaScores <- getScores.SWATH(MS2, input$tarID, msDB, ppm=50, adduct='M+H', typeDB='experimental', eval='median')

