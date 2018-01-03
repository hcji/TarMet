data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

function(input, output){
  options(shiny.maxRequestSize=1024^4)

  # Isotopic Analysis
  raw <- reactive({
    req(input$file)
    withProgress(message = 'Loading Data', value = 0.3, {
      LoadData(input$file$datapath)
    })
  })
  output$rtControl <-  renderUI({
    tagList(
      numericInput('rtleft', 'Input the start retention time', min(raw()$times)),
      numericInput('rtright', 'Input the end retention time', max(raw()$times))
    )
  })
  output$formula_contral <- renderUI({
    if (input$target_select == 'formula') {
      tagList(
        textInput('formula', 'Input the targeted metabolite'),
        selectInput('adduct', 'Select the type of adduct', choices = list(
          Positive = adducts$Name[adducts$Ion_mode == 'positive'],
          Negative = adducts$Name[adducts$Ion_mode == 'negative']
        ))
      )
    }else if (input$target_select == 'm/z of ion'){
      tagList(
        numericInput('fmz', 'Input the monoisotopic mass of the targeted ion of metabolite', 0),
        numericInput('fcharge', 'Input the number of charge of ion', 1)
      )
    }
  })
  
  EICs <- eventReactive(input$button_formula,{
    raw <- raw()
    if (input$target_select == 'formula') {
      eics <- getIsoEIC.formula(raw, input$formula, nmax = input$nmax, adduct = input$adduct, ppm = input$ppm, rtrange = c(input$rtleft, input$rtright), threshold = input$threshold)
    } else if (input$target_select == 'm/z of ion'){
      eics <- getIsoEIC.mz(raw, input$fmz, nmax = input$nmax, ppm = input$ppm, rtrange = c(input$rtleft, input$rtright), charge = input$fcharge)
    }
    
    if (input$ifsmooth){
      eics$eics <- lapply(eics$eics, function(eic){
        if (sum(eic$intensity) > 0) {eic$intensity <- eic$intensity - airPLS(eic$intensity)}
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
    ind <- which.max(peaks$PeakArea[1,-c(1,2)])
    
    tagList(
      numericInput('target_position', 'Define the target retention time: peak position',  peaks$PeakInfo$Position[ind] ),
      numericInput('target_left', 'Define the target retention time: start',  peaks$PeakInfo$Start[ind] ),
      numericInput('target_right', 'Define the target retention time: end',  peaks$PeakInfo$End[ind] )
    )
  })
  
  output$EICPlot <- renderPlotly({
    eics <- EICs()
    peaks <- Peaks()
    withProgress(message = 'Creating plot', value = 0.5, {
      p <- plotEICs(eics, peaks)
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
  
  # Quantitative Analysis
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
    withProgress(message = 'Loading Data', value = 0.3, {
      rawfiles <- lapply(filepathes(), LoadData)
      setProgress(0.6)
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
        intensity.new <- approx(x=eics[[i]]$rt, y=eics[[i]]$intensity, xout=rt, rule = 2)$y
        intensity.new[is.na(intensity.new)] <- 0
        if (input$ifsmooth){
          intensity.new <- intensity.new - airPLS(intensity.new)
        }
        eics[[i]]$rt <- rt
        eics[[i]]$intensity <- intensity.new
      }
      setProgress(0.9)
      if (input$ifalign){
        for (i in 1:length(eics)){
          intensity.new <- PAFFT(t(eics[[i]]$intensity), t(eic$intensity), align.seg, align.shift)$alignedSpectrum
          eics[[i]]$intensity <- as.numeric(intensity.new)
        }
      }
    })
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
    withProgress(message = 'Creating plot', value = 0.5, {
      p <- plotEICs(list(eics=eics), peaks = NULL, Names = Names)
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
  
  # Swath Analysis
  SwathRaw <- reactive({
    withProgress(message = 'Loading Data', value = 0.5, {
      preMz <- mean(EICs()$mzs[1,])
      LoadSwath(input$file$datapath, preMz)
    })
  })
  
  SwathRaws <- reactive({
    req(input$files)
    res <- list()
    preMz <- mean(EICs()$mzs[1,])
    withProgress(message = 'Loading Data', value = 0.1, {
      for (i in seq_along(input$files$datapath)){
        incProgress(1/length(input$files$datapath))
        res[[i]] <- LoadSwath(input$files$datapath[i], preMz)
      }
    })
    res
  })
  
  SwathMz <- reactive({
    raw_swath <- SwathRaw()
    getSwathMz(raw_swath, input$target_position, input$ppm)
  })
  
  SwathEICs <- reactive({
    raw_swath <- SwathRaw()
    swath_mz <- SwathMz()
    getSwathEICs(raw_swath, swath_mz, input$ppm)
  })
  
  SwathCharIon <- reactive({
    swath_eics <- SwathEICs()
    swath_mz <- SwathMz()
    eics <- EICs()
    eic_ms1 <- eics$eics[[1]]
    getSwathCharIon(swath_eics, swath_mz, eic_ms1, input$target_left, input$target_right, input$corr_thres, input$int_thres)
  })
  
  SwathArea <- reactive({
    raw_swath <- SwathRaw()
    raw_swaths <- SwathRaws()
    char_ion <- SwathCharIon()
    rt <- EICs()$eics[[1]]$rt
    align.seg <- round(input$align.seg / mean(diff(rt)))
    align.shift <- round(input$align.shift / mean(diff(rt)))
    withProgress(message = 'Calculating Areas of MS2 EICs', value = 0.5, {
      res <- getSwathArea(raw_swaths, char_ion, input$ppm, input$target_left, input$target_right, align.seg, align.shift, input$area_thres)
    })
    res
  })
  
  output$SwathEicPlot <- renderPlotly({
    swath_area <- SwathArea()
    eics <- swath_area$swaths_eics[[input$swath_ind]]$eics
    p <- plotEICs(list(eics=eics), Names = round(swath_area$swath_mz, 3))
    p <- add_markers(p, x = c(input$target_left, input$target_right), y = rep(0, 2), name = 'peak bound', color = I('blue'), marker = list(size = 5))
    p
  })
  
  output$SwathAreaTable <- renderTable({
    swath_area <- SwathArea()
    Names <- sapply(input$files$name, function(f){
      Name <- strsplit(f,'/')[[1]]
      Name[length(Name)]
    })
    as.data.frame(cbind(Names , round(swath_area$swaths_areas,2)))
  })
  
  output$SwathMS2 <- renderTable({
    swath_area <- SwathArea()$
    rel_area <- swath_area$swaths_areas / swath_area$swaths_areas[,1]
    ms2_int <- round(colMeans(rel_area)/max(colMeans(rel_area)) * 100, 2)
    ms2_mz <- colnames(rel_area)
    data.frame(mz=ms2_mz, abundance=ms2_int)
  })
  
}