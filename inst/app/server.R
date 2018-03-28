library(shiny)
library(plotly)

data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

function(input, output){
  # set space
  options(shiny.maxRequestSize=1024^4)
  
  # input UI
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
  
  output$targetCtrl1 <- renderUI({
    if (input$input=='config file') {
      tagList(
        fileInput('config', 'Choose the config file')
      )
    } else {
      if (input$define=='mass-to-charge'){
        tagList(
          numericInput('mz1', 'Targeted mass-to-charge:', 0)
        ) 
      } else {
        tagList(
          textInput('formula1', 'Targeted compound:')
        ) 
      }
    }
  })
  
  output$targetCtrl2 <- renderUI({
    req(input$config)
    choice <- getTargets(input$config$datapath, input$define)
    if (input$define=='mass-to-charge'){
      tagList(
        selectInput('mz2', 'Targeted mass-to-charge:', choice=choice)
      ) 
    } else {
      tagList(
        selectInput('formula2', 'Targeted compound:', choice=choice)
      ) 
    }
  })
  
  output$formulaCtrl <- renderUI({
    if (input$define=='mass-to-charge'){
      if (input$input=='config file'){
        mz <- as.numeric(input$mz2)
      } else {
        mz <- as.numeric(input$mz1)
      }
      choice <- try(getCandidateFormula(mz, window=0.005, adduct=input$adduct))
      if (class(choice)=='try-error'){
        choice <- NULL
      }
      tagList(
        selectInput('formula3', 'Select the possible formula', choice=choice)
      ) 
    }
  })
  
  formula <- reactive({
    req(input$files)
    if (input$input=='config file' && input$define!='mass-to-charge') {
      as.character(input$formula2)
    } else if (input$input=='direct' && input$define!='mass-to-charge'){
      as.character(input$formula1)
    } else {
      as.character(input$formula3)
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
    getIsotopicPattern(formula(), input$adduct, input$threshold, input$resolution)
  })
  
  targetMzRanges <- reactive({
    if (input$type=='isotopic tracer'){
      mzs <- getMzWithTracer(formula(), input$adduct, input$tracer_element, input$tracer_isotope, input$tracer_number)
    } else {
      mzs <- pattern()[,1]
    }
    getMzRanges(mzs, input$resolution)
  })
  
  # load raw dataset
  rawDataset <- reactive({
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
  
  targetEICs <- eventReactive(input$confirm, {
    res <- list()
    withProgress(message = 'Generating EICs', value = 0.1, {
      for (i in seq_along(sampleNames())){
        res[[i]] <- getMzEICs(rawDataset()[[i]], mzranges=targetMzRanges(), baseline=input$baseline, smooth=input$smooth)
      }
      incProgress(1/length(input$files$datapath))
    })
    names(res) <- sampleNames()
    res
  })
  
  targetPeaks <- reactive({
    req(targetEICs())
    if (input$type!='isotopic tracer'){
      theoretical <- as.numeric(pattern()[,2])
    } else {
      theoretical <- NULL
    }
    getIsotopicPeaks(targetEICs()[[sampleInd()]], SNR.Th=input$snr.th, peakScaleRange=input$scale.th, peakThr=input$int.th, theoretical=theoretical)
  })
  
  output$EICPlot <- renderPlotly({
    plotEICs(targetEICs()[[sampleInd()]], targetPeaks())
  })
  
  whichPeak <- reactive({
    if (input$type!='isotopic tracer'){
      which.max(targetPeaks()$PeakInfo$Similarity)
    } else {
      which.max(colSums(targetPeaks()$PeakArea))
    }
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
      Similarity <- round(getIsotopicSimilarity(userArea(), pattern()[,2]), 3)
      data.frame(Name='User', Position=Position, Start=Start, End=End, Similarity=Similarity)
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
    res <- cbind(res, userArea())
    colnames(res)[1] <- 'mzRange'
    colnames(res)[ncol(res)] <- 'User'
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
      getAlignedEICs(targetEICs(), sampleInd(), input$shift, input$segment)
    } else {
      targetEICs()
    }
  })
  
  output$SamplesPlot <- renderPlotly({
    EicstoPlot <- lapply(alignedEICs(), function(eics){
      eics[[input$WhtoPlot]]
    })
    plotEICs(EicstoPlot, rtrange=c(input$targetRtLeft, input$targetRtRight))
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
      formula <- formula()
      targetMzRanges <- targetMzRanges()
      rtleft <- input$targetRtLeft
      rtright <- input$targetRtRight
      
      quant_res_this <- cbind(formula, targetMzRanges, rtleft, rtright, SamplesArea)
      quant_res(rbind(quant_res(), quant_res_this))
    })
  })
  
  output$resultList <- renderTable({
    quant_res()
  })
  
  output$resultDown <- downloadHandler(
    filename = "results.csv",
    content = function(file) {
      write.csv(quant_res(), filename, row.names = FALSE)
    },
    contentType = "text/csv"
  )

}
