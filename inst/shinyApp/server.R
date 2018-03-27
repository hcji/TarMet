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
      selectInput('sample', 'Select a sampe as reference', choices = sampleNames())
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
          numericInput('mz1', 'Input the targeted mass-to-charge', 0)
        ) 
      } else {
        tagList(
          textInput('formula1', 'Input the targeted compound')
        ) 
      }
    }
  })
  
  output$targetCtrl2 <- renderUI({
    req(input$config)
    choice <- getTargets(input$config$datapath, input$define)
    if (input$define=='mass-to-charge'){
      tagList(
        selectInput('mz2', 'Select the targeted mass-to-charge', choice=choice)
      ) 
    } else {
      tagList(
        selectInput('formula2', 'Select the targeted compound', choice=choice)
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
  
  formula <- eventReactive(input$confirm, {
    if (input$input=='config file' && input$define!='mass-to-charge') {
      as.character(input$formula2)
    } else if (input$input=='direct' && input$define!='mass-to-charge'){
      as.character(input$formula1)
    } else {
      as.character(input$formula3)
    }
  })
  
  output$tracerCtrl1 <- renderUI({
    if (input$type=='isotopic tracer'){
      tagList(
        selectInput('tracer_element', 'Select the element of isotopic tracer', getElements(formula()))
      )
    } else {
      tagList(
        numericInput('threshold', 'Input the threshold of abundance of isotopic peaks', 0.01)
      )
    }
  })
  
  output$tracerCtrl2 <- renderUI({
    req(input$tracer_element)
    tagList(
      selectInput('tracer_isotope', 'Select the isotope of tracer', isotopes$isotope[isotopes$element == input$tracer_element][-1]),
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
  
  targetEICs <- reactive({
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
  
  output$PeakInfo <- renderTable({
    targetPeaks()$PeakInfo
  })
  
  output$PeakArea <- renderTable({
    targetPeaks()$PeakArea
  })
  
  whichPeak <- reactive({
    which.max(targetPeaks()$PeakInfo$Similarity)
  })
  
  output$targetRtCtrl <- renderUI({
    tagList(
      h4('Targeted Retention Time'),
      numericInput('targetRtLeft', 'Targeted retention time (Left)', targetPeaks()$PeakInfo$Start[whichPeak()]),
      numericInput('targetRtRight', 'Targeted retention time (Reft)', targetPeaks()$PeakInfo$End[whichPeak()])
    ) 
  })
  

}
