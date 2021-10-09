library(shiny)
library(plotly)

data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

function(input, output){
  # set space
  options(shiny.maxRequestSize=1024^4)
  
  # load dataset
  rawDataset <- reactive({
    req(input$files)
    res <- list()
    withProgress(message = 'Reading Data', value = 0.1, {
      for (i in seq_along(sampleNames())){
        res[[i]] <- LoadData(input$files$datapath[i])
      }
      incProgress(1/length(input$files$datapath))
    })
    names(res) <- sampleNames()
    res
  })
  
  # define samples
  sampleNames <- reactive({
    req(input$files)
    getSampleName(input$files$name)
  })
  
  output$sampleCtrl <- renderUI({
    req(input$files)
    tagList(
      selectInput('sample', 'Reference Sample', choices = sampleNames())
    )
  })
  
  sampleInd <- reactive({
    which(input$sample==sampleNames())
  })
  
  # define targets
  output$targetCtrl1 <- renderUI({
    if (input$input=='config file') {
      tagList(
        fileInput('config', 'Choose the csv config file', accept = ".csv")
      )
    } else {
      tagList(
        textInput('compoundName', 'Name of compound'),
        selectInput('define', 'How to define the new compound', c('formula', 'mass-to-charge'))
      )
    }
  })
  
  # if define targeted compound via config file
  config <- reactive({
    if (input$input=='config file') {
      req(input$config)
      read.csv(input$config$datapath)
    }
  })
  
  output$targetCtrl2 <- renderUI({
    if (input$input=='config file') {
      req(config())
      config.name <- as.character(config()$name)
      tagList(
        selectInput('target', 'Select a targeted compound', choices = config.name)
      )
    }
  })
  
  default <- reactive({
    config <- config()
    if (input$input=='config file'){
      wh <- which(config$name==input$target)
      formula <- if (!is.na(config$formula[wh])) {config$formula[wh]} else {''}
      adduct <- if (!is.na(config$adduct[wh])) {config$adduct[wh]} else {'M+H'}
      ppm <- if (!is.na(config$ppm[wh])) {config$ppm[wh]} else {50}
      rtmin <- if (!is.na(config$rtmin[wh])) {config$rtmin[wh]} else {0}
      rtmax <- if (!is.na(config$rtmax[wh])) {config$rtmax[wh]} else {max(rawDataset()[[1]]$times)}
      scale <- if (!is.na(config$scale[wh])) {config$scale[wh]} else {5}
      height <- if (!is.na(config$height[wh])) {config$height[wh]} else {0}
      snr <- if (!is.na(config$snr[wh])) {config$snr[wh]} else {5}
    } else {
      formula <- ''
      adduct <- 'M+H'
      ppm <- 10
      rtmin <- 0
      rtmax <- max(rawDataset()[[1]]$times)
      scale <- 5
      height <- 0
      snr <- 5
    }
    list(formula=formula, adduct=adduct, ppm=ppm, rtmin=rtmin, rtmax=rtmax, scale=scale, height=height, snr=snr)
  })
  
  output$paraCtrl0 <- renderUI({
    if (input$input != 'config file'){
      if (input$define == 'formula'){
        tagList(
          textInput('formula', 'Formula of target:', value = default()$formula),
          selectInput('adduct', 'Adduct type:', choices = list(
            Positive = adducts$Name[adducts$Ion_mode == 'positive'],
            Negative = adducts$Name[adducts$Ion_mode == 'negative']
          ), selected = default()$adduct)
        )
      }else{
        tagList(
          numericInput('mass_to_charge', 'mass-to-charge:', value = 0),
          numericInput('charge', 'mass-to-charge:', value = 1)
        )
      }
    }else{
      tagList(
        textInput('formula', 'Formula of target:', value = default()$formula),
        selectInput('adduct', 'Adduct type:', choices = list(
          Positive = adducts$Name[adducts$Ion_mode == 'positive'],
          Negative = adducts$Name[adducts$Ion_mode == 'negative']
        ), selected = default()$adduct)
      )
    }
  })
  
  output$paraCtrl <- renderUI({
    tagList(
      h4('EIC extraction'),
      numericInput('ppm', 'EIC ppm:', default()$ppm),
      numericInput('rtmin', 'RT Start:', default()$rtmin),
      numericInput('rtmax', 'RT End:', default()$rtmax),
      h4('Peak Detection'),
      selectInput('baseline', 'Remove Baseline?', c(FALSE, TRUE)),
      selectInput('smooth', 'Smooth EIC?', c(FALSE, TRUE)),
      selectInput('fineness', 'Fineness of Detection', c('Medium', 'High', 'Low')),
      numericInput('snr.th', 'SNR threshold', default()$snr),
      numericInput('scale.th', 'Scale threshold (s)', default()$scale),
      numericInput('int.th', 'Intensity threshold (above the baseline)', default()$height)
    )
  })
  
  # define formula
  formula <- reactive({
    gsub(" ", "", as.character(input$formula))
  })
  
  compoundName <- reactive({
    if (input$input=='config file') {
      input$target
    } else {
      input$compoundName
    }
  })
  
  output$tracerCtrl1 <- renderUI({
    req(formula())
    if (input$type=='isotopic tracer'){
      tagList(
        selectInput('tracer_element', 'Element of isotopic tracer:', getElements(formula()))
      )
    } else {
      tagList(
        numericInput('threshold', 'Threshold of abundance of isotopic peaks', 0.01)
      )
    }
  })
  
  output$tracerCtrl2 <- renderUI({
    req(input$tracer_element)
    tagList(
      selectInput('tracer_isotope', 'Isotope of tracer', isotopes$isotope[isotopes$element == input$tracer_element][-1]),
      numericInput('tracer_number', 'Input n, where at most M+n isotopologues are detected.', min(3, getElementNum(formula(), input$tracer_element)), 
                   min = 1, max = getElementNum(formula(), input$tracer_element))
    )
  })
  
  # generate target compound informations
  pattern <- reactive({
    req(formula())
    isolate(getIsotopicPattern(formula(), input$adduct, input$threshold, input$resolution))
  })
  
  targetMzRanges <- reactive({
    if (input$type=='isotopic tracer'){
      mzs <- getMzWithTracer(formula(), input$adduct, input$tracer_element, input$tracer_isotope, input$tracer_number)
    } else if(input$define == 'formula'){
      mzs <- pattern()[,1]
    } else {
      mzs <- c(input$mass_to_charge) + 1.0034 * c(0, 1, 2, 3) / input$charge
    }
    isolate(getMzRanges(mzs, resolution=input$resolution, ppm=input$ppm))
  })
  
  targetEICs <- eventReactive(input$confirm, {
    res <- list()
    rtranges <- c(input$rtmin, input$rtmax)
    isolate(withProgress(message = 'Generating EICs', value = 0.1, {
      for (i in seq_along(sampleNames())){
        res[[i]] <- getMzEICs(rawDataset()[[i]], rtranges=rtranges, mzranges=targetMzRanges(), baseline=input$baseline, smooth=input$smooth)
      }
      incProgress(1/length(input$files$datapath))
    }))
    names(res) <- sampleNames()
    res
  })
  
  targetPeaks <- reactive({
    req(targetEICs())
    isolate(if (input$type!='isotopic tracer'){
      if (input$define == 'formula'){
        theoretical <- as.numeric(pattern()[,2])
      } else {
        theoretical <- NULL
      }
    } else {
      theoretical <- NULL
    })
    withProgress(message = 'Peak detection', value = 0.5, {
      getIsotopicPeaks(targetEICs()[[sampleInd()]], SNR.Th=input$snr.th, peakScaleRange=input$scale.th, peakThr=input$int.th, theoretical=theoretical, fineness=input$fineness)
    })
  })
  
  output$EICPlot <- renderPlotly({
    plotEICs(targetEICs()[[sampleInd()]], targetPeaks())
  })
  
  whichPeak <- reactive({
    isolate(if (input$type!='isotopic tracer'){
      if (input$define == 'formula'){
        which.max(targetPeaks()$PeakInfo$Similarity)
      } else {
        which.max(colSums(targetPeaks()$PeakArea))
      }
    } else {
      which.max(colSums(targetPeaks()$PeakArea))
    })
  })
  
  output$targetRtCtrl <- renderUI({
    tagList(
      h4('Targeted Retention Time'),
      numericInput('targetRtLeft', 'Targeted retention time (Left)', targetPeaks()$PeakInfo$Start[whichPeak()]),
      numericInput('targetRtRight', 'Targeted retention time (Reft)', targetPeaks()$PeakInfo$End[whichPeak()]),
      downloadButton("targetDown", "Download")
    ) 
  })
  
  userArea <- reactive({
    round(getArea(targetEICs()[[sampleInd()]], input$targetRtLeft, input$targetRtRight),2)
  })
  
  userInfo <- reactive({
    Position <- NA
    Start <- input$targetRtLeft
    End <- input$targetRtRight
    if (input$type!='isotopic tracer'){
      if (input$define == 'formula'){
        Similarity <- round(getIsotopicSimilarity(userArea(), pattern()[,2]), 3)
        data.frame(Name='User', Position=Position, Start=Start, End=End, Similarity=Similarity)
      } else {
        data.frame(Name='User', Position=Position, Start=Start, End=End)
      }
    } else {
      data.frame(Name='User', Position=Position, Start=Start, End=End)
    }
  })
  
  outputPeakInfo <- reactive({
    rbind(targetPeaks()$PeakInfo[,1:ncol(userInfo())], userInfo())
  })
  
  outputPeakArea <- reactive({
    res <- targetPeaks()$PeakArea
    res <- cbind(rownames(res), res)
    res <- cbind(res, userArea(), round(userArea()/sum(userArea()),3))
    colnames(res)[1] <- 'mzRange'
    colnames(res)[(ncol(res)-1) : ncol(res)] <- c('User','Abundance')
    res
  })
  
  output$PeakInfo <- renderTable(outputPeakInfo())
  output$PeakArea <- renderTable(outputPeakArea())
  
  output$targetDown <- downloadHandler(
    filename = paste(formula(), 'csv', sep='.'),
    content = function(filename) {
      write.csv(outputPeakArea(), filename, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Alignment
  output$alignmentCtrl <- renderUI({
    if (input$alignment){
      tagList(
        numericInput('shift', 'Maximum retention time of shift (s)', 20),
        numericInput('segment', 'Size of segment of PAFFT (s)', 20)
      )
    }
  })
  
  output$isoPlotCtrl <- renderUI({
    tagList(
      selectInput('WhtoPlot', 'Which isotopic peak to plot', choices = as.character(rownames(targetPeaks()$PeakArea)))
    )
  })
  
  alignedEICs <- reactive({
    if (input$alignment){
      getAlignedEICs(targetEICs(), sampleInd(), input$shift, input$segment, align=TRUE)
    } else {
      getAlignedEICs(targetEICs(), sampleInd(), input$shift, input$segment, align=FALSE)
    }
  })
  
  output$SamplesPlot <- renderPlotly({
    withProgress(message = 'Generating Plots', value = 0.5, {
      EicstoPlot <- lapply(alignedEICs(), function(eics){
        eics[[input$WhtoPlot]]
      })
      plotEICs(EicstoPlot, rtrange=c(input$targetRtLeft, input$targetRtRight))
    })
  })
  
  outputSamplesArea <- reactive({
    getQuantifiedResult(alignedEICs(), input$targetRtLeft, input$targetRtRight)
  })
  
  output$SamplesArea <- renderTable({
    res <- outputSamplesArea()
    res <- cbind(rownames(res), res)
    colnames(res)[1] <- 'mzRange'
    res
  })
  
  output$SamplesBarPlot <- renderPlotly({
    plotStackBar(outputSamplesArea())
  })
  
  # Result List
  quant_res <- reactiveVal(data.frame())
  
  observeEvent(input$addButton, {
    withProgress(message = 'Successful', value = 0.3, {
      SamplesArea <- outputSamplesArea()
      compound <- compoundName()
      formula <- formula()
      targetMzRanges <- targetMzRanges()
      rtleft <- input$targetRtLeft
      rtright <- input$targetRtRight
      
      quant_res_this <- cbind(compound, formula, targetMzRanges, rtleft, rtright, SamplesArea)
      quant_res(rbind(quant_res(), quant_res_this))
    })
  })
  
  output$resultList <- renderTable({
    quant_res()
  })
  
  output$resultDown <- downloadHandler(
    filename = "results.csv",
    content = function(filename) {
      write.csv(quant_res(), filename, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  config.new <- reactiveVal(data.frame())
  observeEvent(input$addConfig, {
    name <- compoundName()
    formula <- formula()
    mz <- pattern()[1,1]
    adduct <- input$adduct
    ppm <- input$ppm
    rtmin <- input$rtmin
    rtmax <- input$rtmax
    scale <- input$scale.th
    height <- input$int.th
    snr <- input$snr.th
    
    this.config <- cbind(name, formula, mz, adduct, ppm, rtmin, rtmax, scale, height, snr)
    config.new(rbind(config.new(), this.config))
  })
  
  output$configList <- renderTable({
    config.new()
  })
  
  output$configDown <- downloadHandler(
    filename = "config.csv",
    content = function(filename) {
      write.csv(config.new(), filename, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # MS/MS processing
  compMSMS <- reactive({
    i <- sampleInd()
    rawDataMS2 <- LoadMSMS(input$files$datapath[i])
    precursorMzRange <- as.numeric(targetMzRanges()[1,])
    precursorRtRange <- c(input$targetRtLeft, input$targetRtRight)
    getMSMS(rawDataMS2$header, rawDataMS2$peaks, precursorMzRange, precursorRtRange, ppm = input$MatchPPM)
  })
  
  output$compMSMSPlot <- renderPlotly({
    withProgress(message = 'Reading Data', value = 0.5, {
      ms2 <- compMSMS()
      if (!is.null(ms2)){
        plotMS(ms2)
      }
    })
  })
  
  MSMSInfo <- reactive({
    info <- c(compoundName(),
              formula(),
              input$adduct,
              mean(as.numeric(targetMzRanges()[1,])),
              mean(c(input$targetRtLeft, input$targetRtRight))
    )
    names(info) <- c('NAME', 'FORMULA', 'PRECURSORTYPE', 'PRECURSORMZ', 'RETENTIONTIME')
    info
  })
  
  output$MSMSDown <- downloadHandler(
    filename = paste(compoundName(), ".msp", sep=''),
    content = function(filename) {
      WriteMSP(compID, list(info = MSMSInfo(), spec=compMSMS()), filename)
    },
    contentType = "msp"
  )
  
}



