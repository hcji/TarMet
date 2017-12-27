data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

function(input, output, session){
  options(shiny.maxRequestSize=1024^4)
  
  session$onSessionEnded(function() {
    stopApp()
    q("no")
  })
  
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
  
  EICs <- eventReactive(input$button_formula,{
    raw <- raw()
    eics <- getIsoEIC(raw, input$formula, input$fmz, input$nmax, adduct = input$adduct, ppm = input$ppm, rtrange = c(input$rtleft, input$rtright), threshold = input$threshold)
    if (input$ifsmooth){
      eics$eics <- lapply(eics$eics, function(eic){
        if (sum(eic$intensity) > 0) {eic$intensity <- eic$intensity - airPLS(eic$intensity)}
        eic
      })
    }
    eics
  })
  
  EIC <- reactive({
    eics <- EICs()
    ind <- which.max(sapply(eics$eics, function(s){
      sum(s$intensity)
    }))
    eics$eics[[ind]]
  })
  
  Peaks <- reactive({
    eics <- EICs()
    getIsoPeaks(eics, input$SNR.Th, input$peakScaleRange, input$peakThr)
  })
  
  output$targetControl <-  renderUI({
    peaks <- Peaks()
    eics <- EICs()
    rt <- eics$eics[[1]]$rt
    ind <- which.max(peaks$PeakArea[1,-c(1,2)])
    
    tagList(
      numericInput('target_left', 'Define the target retention time: start',  peaks$PeakInfo$Start[ind] ),
      numericInput('target_right', 'Define the target retention time: end',  peaks$PeakInfo$End[ind] )
    )
  })
  
  output$EICPlot <- renderPlotly({
    eics <- EICs()
    peaks <- Peaks()
    eic <- EIC()
    withProgress(message = 'Creating plot', value = 0.1, {
      p <- plot_ly() %>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'Intensity'))
      incProgress(0.1)
      for (f in 1: length(eics$eics)) {
        p <- add_trace(p, x = eics$eics[[f]]$rt, y = eics$eics[[f]]$intensity, mode='line', name = paste('Mz: ',round(eics$mzs[f,1],4), ' - ', round(eics$mzs[f,2],4)))
        incProgress(0.1)
      }
      if (!is.null(peaks)){
        p <- add_markers(p, x = eic$rt[peaks$Index$Position], y = eic$intensity[peaks$Index$Position], name = 'peak position', color = I('red'), marker = list(size = 5))
        p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = eic$intensity[c(peaks$Index$Start, peaks$Index$End)], name = 'peak bound', color = I('blue'), marker = list(size = 5))
      }
    })
    p
  })
  
  output$PeakArea <- renderTable({
    peaks <- Peaks()
    eics <- EICs()
    UserPeakArea <- getArea(eics, input$target_left, input$target_right)
    Ratio <- UserPeakArea/UserPeakArea[1] * 100
    PeakArea <- cbind(peaks$PeakArea, UserPeakArea, Ratio)
    colnames(PeakArea)[((ncol(PeakArea)-1) : ncol(PeakArea))] <- c('User Define', 'relative area (user)')
    
    if(input$fmz < 0){
      pattern <- getIsoPat(input$formula, input$adduct, threshold = input$threshold)
      nmax <- min(input$nmax,round(max(pattern[,1] - pattern[1,1])))
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
  
  output$EICs_Control <-  renderUI({
    mzrange <- EICs()$mzs
    selection <- (1:nrow(mzrange))
    tagList(
      selectInput('eics_iso','Select which isotope feature is used for quantification',selection)
    )
  })
  
  files_eics <- reactive({
    rawfiles <- lapply(filepathes(), LoadData)
    ind <- as.numeric(input$eics_iso)
    eic <- EICs()$eics[[ind]]
    mzrange <- EICs()$mzs[ind,]
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
    Names <- sapply(input$files$name, function(f){
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
        p <- add_trace(p, x = eics[[f]]$rt, y = eics[[f]]$intensity, mode='line', name = paste(Names[f]))
        incProgress(0.1)
      }
      p <- add_markers(p, x = eic$rt[c(peaks$Index$Start, peaks$Index$End)], y = rep(0, 2*length(peaks$Index$End)), name = 'peak bound', color = I('blue'), marker = list(size = 5))
    })
    p
  })
  
  output$files_peaks <- renderTable({
    files_eics <- files_eics()
    rt <- EICs()$eics[[1]]$rt
    Names <- sapply(input$files$name, function(f){
      Name <- strsplit(f,'/')[[1]]
      Name[length(Name)]
    })
    Areas <- sapply(files_eics, function(eics){
      eics <- list(rt=rt, eics=list(eics))
      getArea(eics, input$target_left, input$target_right)
    })
    Areas <- round(Areas, 3)
    res <- cbind(Names, Areas)
    colnames(res) <- c('Names', paste('Area (M+', as.numeric(input$eics_iso)-1, ')', sep='')) 
    res
  })
}