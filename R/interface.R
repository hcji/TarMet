runTarMet <- function(){
  library(shiny)
  library(plotly)
  library(enviPat)
  library(MassSpecWavelet)
  library(baselineWavelet)
  library(data.table)
  library(Matrix)
  
  data("isotopes", package = "enviPat")
  data("adducts", package = "enviPat")
  
  ui <- navbarPage("Targeted Bio-Mass Analyzer",
    tabPanel(
      'Isotopic Analysis',
      titlePanel("Isotopic Analysis"),
      sidebarLayout(
        sidebarPanel(
          fileInput('file', 'Choose one of your raw data files'),
          
          h4('Metabolite Information'),
          textInput('formula', 'Input the targeted metabolite'),
          numericInput('fmz', 'Input the m/z of the targeted metabolite (if the formula is unknow)', -1),
          selectInput('adduct', 'Select the type of adduct', choices = list(
            Positive = adducts$Name[adducts$Ion_mode == 'positive'],
            Negative = adducts$Name[adducts$Ion_mode == 'negative']
          )),
          
          h4('Isotopic Information (only used when the formula is given)'),
          numericInput('threshold', 'Input the threshold of the relative abundance of isotopic peaks', 0.001),
          
          h4('EIC Information'),
          numericInput('ppm', 'Input the tolerance of the m/z difference (ppm)', 100),
          uiOutput('rtControl'),
          selectInput('ifsmooth', 'Select wether the EICs are smoothed or not', c(TRUE, FALSE)),
          
          h4('Peak Detection Information'),
          numericInput('SNR.Th', 'Input the minimum SNR of peaks', 5),
          numericInput('peakScaleRange', 'Input the scale range of the peak', 5),
          numericInput('peakThr', 'Input the minimal absolute intensity (above the baseline) of peaks to be picked', 0),
          uiOutput('targetControl')
          
        ),
        
        mainPanel(
          h4('Extracted EICs'),
          plotlyOutput('EICPlot'),
          h4('Peak Information'),
          tableOutput('PeakInfo'),
          tableOutput('PeakArea')
        )
      )
    ),
    
    tabPanel(
      "Quantitative Analysis",
      titlePanel("Quantitative Analysis"),
      sidebarLayout(
        sidebarPanel(
          fileInput('files', 'Choose multiple raw data files', multiple=TRUE),
          
          h4('Alignment Information'),
          selectInput('ifalign', 'Select wether the EICs are aligned across samples or not', c(TRUE, FALSE)),
          numericInput('align.shift', 'Input the maximum scans of shift', 20),
          numericInput('align.seg', 'Input the size of segment of PAFFT', 20),
          
          h5('The other parameters are the same as the isotopic analysis step. If you want to change them, please go back to the last step.')
        ),
        
        mainPanel(
          h4('Extracted EICs'),
          plotlyOutput('EICPlots'),
          h4('Quantitative Results'),
          tableOutput('files_peaks')
        )
      )
      
    )
  )
  
  server <- function(input, output){
    options(shiny.maxRequestSize=1024^3)
    
    raw <- reactive({
      req(input$file)
      LoadData(input$file$datapath)
    })
    output$rtControl <-  renderUI({
      tagList(
        numericInput('rtleft', 'Input the start retention time', min(raw()$times)),
        numericInput('rtright', 'Input the end retention time', max(raw()$times))
      )
    })

    EICs <- reactive({
      raw <- raw()
      adduct <- which(input$adduct == adducts$Name)
      eics <- getIsoEIC(raw, input$formula, input$fmz, adduct = adduct, ppm = input$ppm, rtrange = c(input$rtleft, input$rtright), threshold = input$threshold)
      if (input$ifsmooth){
        eics$eics <- lapply(eics$eics, function(eic){
          eic$intensity <- eic$intensity - airPLS(eic$intensity)
          eic
        })
      }
      eics
    })
    
    Peaks <- reactive({
      eics <- EICs()
      getIsoPeaks(eics, input$SNR.Th, input$peakScaleRange, input$peakThr)
    })
    
    output$targetControl <-  renderUI({
      peaks <- Peaks()
      eics <- EICs()
      rt <- eics$eics[[1]]$rt
      ind <- which.max(peaks$PeakArea$`Peak 1`)
      
      tagList(
        numericInput('target_left', 'Define the target retention time: start',  peaks$PeakInfo$Start[ind] ),
        numericInput('target_right', 'Define the target retention time: end',  peaks$PeakInfo$End[ind] )
      )
    })
      
    output$EICPlot <- renderPlotly({
      eics <- EICs()
      peaks <- Peaks()
      eic <- eics$eics[[1]]
      withProgress(message = 'Creating plot', value = 0.1, {
        p <- plot_ly() %>%
          layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
                 yaxis = list(title = 'Intensity'))
        incProgress(0.1)
        for (f in 1: length(eics$eics)) {
          p <- add_trace(p, x = eics$eics[[f]]$rt, y = eics$eics[[f]]$intensity, mode='line', name = paste('Mz: ',round(eics$mzs[f,1],4), ' - ', round(eics$mzs[f,2],4)))
          incProgress(0.1)
        }
        p <- add_markers(p, x = eic$rt[peaks$Index$Position], y = eic$intensity[peaks$Index$Position], name = 'peak position', color = I('red'), marker = list(size = 5))
        p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = eic$intensity[c(peaks$Index$Start, peaks$Index$End)], name = 'peak bound', color = I('blue'), marker = list(size = 5))
      })
      p
    })
    
    output$PeakArea <- renderTable({
      peaks <- Peaks()
      eics <- EICs()
      UserPeakArea <- getArea(eics, input$target_left, input$target_right)
      Ratio <- UserPeakArea/UserPeakArea[1]
      PeakArea <- cbind(peaks$PeakArea, UserPeakArea, Ratio)
      colnames(PeakArea)[((ncol(PeakArea)-1) : ncol(PeakArea))] <- c('User Define', 'relative area (user)')
      
      if(input$fmz < 0){
        adduct <- which(input$adduct == adducts$Name)
        pattern <- getIsoPat(formula, adduct, threshold = input$threshold)
        nmax <- round(max(pattern[,1] - pattern[1,1]))
        ints <- sapply(0:nmax, function(n){
          pai <- pattern[round(pattern[,1]- pattern[1,1])==n,]
          sum(pai[,2])
        })
        PeakArea <- cbind(PeakArea, ints)
        colnames(PeakArea)[ncol(PeakArea)] <- 'relative area (theoretical)'
      }
      
      PeakArea
    })
    
    output$PeakInfo <- renderTable({
      peaks <- Peaks()
      peaks$PeakInfo
    })
    
    filepathes <- reactive({
      req(input$files)
      input$files$datapath
    })
    
    files_eics <- reactive({
      eic <- EICs()$eics[[1]]
      rawfiles <- lapply(filepathes(), LoadData)
      mzrange <- EICs()$mzs[1,]
      eics <- lapply(rawfiles, function(raw){
        getEIC(raw, c(input$rtleft, input$rtright), mzrange)
      })
      rt <- EICs()$eics[[1]]$rt
      align.seg <- round(input$align.seg / mean(diff(rt)))
      align.shift <- round(input$align.shift / mean(diff(rt)))
      for (i in 1:length(eics)){
        intensity.new <- approx(x=eics[[i]]$rt, y=eics[[i]]$intensity, xout=rt)$y
        intensity.new[is.na(intensity.new)] <- 0
        if (input$ifsmooth){
          intensity.new <- intensity.new - airPLS(intensity.new)
        }
        eics[[i]]$rt <- rt
        eics[[i]]$intensity <- intensity.new
      }
      if (input$ifalign){
        for (i in 1:length(eics)){
          intensity.new <- PAFFT(t(eics[[i]]$intensity), t(eic$intensity), align.seg, align.shift)$alignedSpectrum
          eics[[i]]$intensity <- as.numeric(intensity.new)
        }
      }
      eics
    })
    
    output$EICPlots <- renderPlotly({
      eics <- files_eics()
      peaks <- Peaks()
      Names <- sapply(filepathes(), function(f){
        Name <- strsplit(f,'/')[[1]]
        Name[length(Name)]
      })
      names(Names) <- NULL
        
      eic <- EICs()$eics[[1]]
      withProgress(message = 'Creating plot', value = 0.1, {
        p <- plot_ly() %>%
          layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
                 yaxis = list(title = 'Intensity'))
        incProgress(0.1)
        for (f in 1: length(eics)) {
          p <- add_trace(p, x = eics[[f]]$rt, y = eics[[f]]$intensity, mode='line', name = paste('Sample: ',Names[f]))
          incProgress(0.1)
        }
        p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = eic$intensity[c(peaks$Index$Start, peaks$Index$End)], name = 'peak bound', color = I('blue'), marker = list(size = 5))
      })
      p
    })
    
    output$files_peaks <- renderTable({
      files_eics <- files_eics()
      rt <- EICs()$eics[[1]]$rt
      Names <- sapply(filepathes(), function(f){
        Name <- strsplit(f,'/')[[1]]
        Name[length(Name)]
      })
      Areas <- sapply(files_eics, function(eics){
        eics <- list(rt=rt, eics=list(eics))
        getArea(eics, input$target_left, input$target_right)
      })
      data.frame(Names, Areas)
    })
    
  }
  
  shinyApp(ui, server)
  
}



