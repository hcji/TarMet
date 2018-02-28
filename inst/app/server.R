data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

function(input, output) {
  
  options(shiny.maxRequestSize=1024^4)
  
  # Isotopologues Analysis
  iso_name_sample <- reactive({
    Names <- sapply(input$iso_files$name, function(f){
      Name <- strsplit(f,'/')[[1]]
      Name[length(Name)]
    })
    as.character(Names)
  })

  output$iso_ctrl_sample <- renderUI({
    req(input$iso_files)
    tagList(
      selectInput('iso_name_sample', 'Select a sampe as reference', choices = iso_name_sample())
    )
  })

  output$iso_ctrl_target <- renderUI({
    if (input$iso_target_select == 'YES'){
      tagList(
        textInput('iso_formula', 'Input the targeted metabolite'),
        selectInput('iso_adduct', 'Select the type of adduct', choices = list(
          Positive = adducts$Name[adducts$Ion_mode == 'positive'],
          Negative = adducts$Name[adducts$Ion_mode == 'negative']
        ))
      )
    } else if (input$iso_target_select == 'NO'){
      tagList(
        numericInput('iso_fmz', 'Input the monoisotopic mass of the targeted ion of metabolite', 0),
        numericInput('iso_fcharge', 'Input the number of charge of ion', 1),
        numericInput('iso_nmax', 'Input n, where at most M+n isotopologues are detected.', 3),
        numericInput('iso_C13', 'Input the maximum number of C13 to be considered.', 3),
        numericInput('iso_H2', 'Input the maximum number of H2 to be considered.', 0),
        numericInput('iso_O18', 'Input the maximum number of O18 to be considered.', 0),
        numericInput('iso_N15', 'Input the maximum number of N15 to be considered.', 0),
        numericInput('iso_S34', 'Input the maximum number of S34 to be considered.', 0)
      )
    }
  })
  
  output$iso_ctrl_tracer1 <- renderUI({
    req(input$iso_formula)
    if (input$iso_target_select == 'YES' && input$iso_assay_type == 'targeted analysis'){
      tagList(
        numericInput('iso_resolution', 'Input the resolution of your MS', 50000),
        numericInput('iso_threshold', 'Input the threshold of the relative abundance of isotopic peaks', 0.001)
      )
    } else if (input$iso_target_select == 'YES' && input$iso_assay_type == 'isotopic tracer'){
      tagList(
        selectInput('iso_tracer_element', 'Select the element of isotopic tracer', getElements(input$iso_formula))
      )
    }
  })
  output$iso_ctrl_tracer2 <- renderUI({
    req(input$iso_tracer_element)
    tagList(
      selectInput('iso_tracer_isotope', 'Select the isotope of tracer', isotopes$isotope[isotopes$element == input$iso_tracer_element][-1]),
      numericInput('iso_tracer_number', 'Input n, where at most M+n isotopologues are detected.', min(3, getElementNum(input$iso_formula, input$iso_tracer_element)), 
                   min = 1, max = getElementNum(input$iso_formula, input$iso_tracer_element))
    )
  })

  iso_raws <- reactive({
    res <- list()
    withProgress(message = 'Reading Data', value = 0.1, {
      for (i in seq_along(iso_name_sample())){
        res[[i]] <- LoadData(input$iso_files$datapath[i])
      }
      incProgress(1/length(input$iso_files$datapath))
    })
    res
  })

  output$iso_ctrl_rt <-  renderUI({
    req(input$iso_files)
    raw_data <- iso_raws()[[1]]
    tagList(
      numericInput('iso_eic_left', 'Input the start retention time', min(raw_data$times)),
      numericInput('iso_eic_right', 'Input the end retention time', max(raw_data$times))
    )
  })

  iso_ref_eic <- eventReactive(input$iso_target_button,{
    withProgress(message = 'Extracting EICs', value = 0.5, {
      ind <- which(iso_name_sample() == input$iso_name_sample)
      raw_data <- iso_raws()[[ind]]
      if (input$iso_target_select == 'YES' && input$iso_assay_type == 'targeted analysis') {
        eics <- getIsoEIC.formula(raw_data, input$iso_formula, adduct = input$iso_adduct,
                                  ppm = input$iso_ppm, rtrange = c(input$iso_eic_left, input$iso_eic_right),
                                  threshold = input$iso_threshold, resolution = input$iso_resolution)
      } else if (input$iso_target_select == 'YES' && input$iso_assay_type == 'isotopic tracer'){
        eics <- getIsoEIC.tracer(raw_data, input$iso_formula, adduct = input$iso_adduct,
                                 ppm = input$iso_ppm, rtrange = c(input$iso_eic_left, input$iso_eic_right),
                                 element = input$iso_tracer_element, isotope = input$iso_tracer_isotope, 
                                 number = input$iso_tracer_number)
      } else if (input$iso_target_select == 'NO'){
        eics <- getIsoEIC.mz(raw_data, input$iso_fmz, nmax = input$iso_nmax, ppm = input$iso_ppm,
                             rtrange = c(input$iso_eic_left, input$iso_eic_right), charge = input$iso_fcharge,
                             elements = c('C','H','O','N','S'), isotope = c('13C', '2H', '18O', '15N', '34S'), 
                             numbers = c(input$iso_C13, input$iso_H2, input$iso_O18, input$iso_N15, input$iso_S34))
      }
      if (input$iso_baseline){
        eics$eics <- lapply(eics$eics, function(eic){
          if (sum(eic$intensity) > 0) {eic$intensity <- eic$intensity - airPLS(eic$intensity)}
          eic
        })
      }
    })
    
    eics
  })

  iso_peaks <- reactive({
    getIsoPeaks(iso_ref_eic(), input$iso_peak_snr.th, input$iso_peak_scale.th, input$iso_peak_int.th)
  })

  output$iso_ctrl_targetRT <-  renderUI({
    peaks <- iso_peaks()
    rt <- iso_ref_eic()$eics[[1]]$rt
    ind <- as.numeric(which.max(peaks$PeakArea[1,-c(1,2)]))

    tagList(
      numericInput('iso_target_position', 'Define the target retention time: peak position', peaks$PeakInfo$Position[ind]),
      numericInput('iso_target_left', 'Define the target retention time: start',  peaks$PeakInfo$Start[ind]),
      numericInput('iso_target_right', 'Define the target retention time: end',  peaks$PeakInfo$End[ind]),
      numericInput('iso_target_isotope', 'Define the target isotopologue index', 1, min = 1, step = 1)
    )
  })

  output$iso_eic_plot <- renderPlotly({
    eics <- iso_ref_eic()
    peaks <- iso_peaks()
    Names <- apply(eics$mzs, 1, function(mz){paste(round(mz[1],4), '-',round(mz[2], 4))})
    withProgress(message = 'Creating EIC Plot', value = 0.5, {
      p <- plotEICs(eics, peaks = peaks, Names = Names)
    })
    p
  })

  iso_peak_area_info <- reactive({
    peaks <- iso_peaks()
    eics <- iso_ref_eic()
    target_area <- getArea(eics, input$iso_target_left, input$iso_target_right)
    target_ratio <- target_area/max(target_area) * 100
    areas <- cbind(peaks$PeakArea, target_area, target_ratio)
    simis <- NULL
    if (input$iso_target_select == 'YES' && input$iso_assay_type == 'targeted analysis') {
      theoretical_ratio <- eics$pattern[,2]
      areas <- cbind(areas, theoretical_ratio)
      
      simis <- sapply(3:ncol(peaks$PeakArea), function(s){
        getIsoSimi(peaks$PeakArea[,s], eics$pattern[,2])
      })
    }
    return(list(areas = areas, simis=simis))
  })
  
  output$iso_peak_info <- renderTable({
    peaks <- iso_peaks()
    if (input$iso_target_select == 'YES' && input$iso_assay_type == 'targeted analysis'){
      peaks$PeakInfo <- cbind(peaks$PeakInfo, iso_peak_area_info()$simis)
      colnames(peaks$PeakInfo)[5] <- 'isotopic_similarity'
    }
    peaks$PeakInfo
  })
  
  output$iso_peak_area <- renderTable({
    iso_peak_area_info()$areas
  })
  
  output$iso_download <- downloadHandler(
    filename = "isotopologues.csv",
    content = function(file) {
      write.csv(iso_peak_area_info()$areas, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  iso_files_eics <- reactive({
    rawfiles <- iso_raws()
    withProgress(message = 'Extracting EICs from files', value = 0.3, {
      ind <- input$iso_target_isotope
      ref_eic <- iso_ref_eic()$eics[[ind]]
      mzrange <- iso_ref_eic()$mzs[ind,]
      rt <- ref_eic$rt
      align.seg <- round(input$iso_align.seg / mean(diff(rt)))
      align.shift <- round(input$iso_align.shift / mean(diff(rt)))
      
      eics <- lapply(rawfiles, function(raw){
        getEIC(raw, c(input$iso_eic_left, input$iso_eic_right), mzrange)
      })
      setProgress(0.6)
      for (i in 1:length(eics)){
        intensity.new <- approx(x=eics[[i]]$rt, y=eics[[i]]$intensity, xout=rt, rule = 2)$y
        if (input$iso_baseline){
          intensity.new <- intensity.new - airPLS(intensity.new)
        }
        eics[[i]]$rt <- rt
        eics[[i]]$intensity <- intensity.new
      }
      setProgress(0.9)
      if (input$iso_align){
        for (i in 1:length(eics)){
          intensity.new <- PAFFT(t(eics[[i]]$intensity), t(ref_eic$intensity), align.seg, align.shift)$alignedSpectrum
          eics[[i]]$intensity <- as.numeric(intensity.new)
        }
      }
    })
    eics
  })

  output$iso_files_EIC_Plot <- renderPlotly({
    files_eics <- iso_files_eics()
    peaks <- iso_peaks()
    Names <- iso_name_sample()
    
    rt <- iso_ref_eic()$eics[[1]]$rt
    withProgress(message = 'Creating plot', value = 0.5, {
      p <- plotEICs(list(eics = files_eics), peaks = NULL, Names = Names)
      p <- add_markers(p, x = rt[c(peaks$Index$Start, peaks$Index$End)], y = rep(0, 2*length(peaks$Index$End)), name = 'peak bound', color = I('blue'), marker = list(size = 5))
    })
    p
  })
  
  iso_files_areas <- reactive({
    files_eics <- iso_files_eics()
    Areas <- sapply(files_eics, function(eics){
      eics <- list(eics=list(eics))
      getArea(eics, input$iso_target_left, input$iso_target_right)
    })
    Areas <- round(Areas, 3)
    Areas
  })
  
  output$iso_files_peaks <- renderTable({
    Names <- iso_name_sample()
    Areas <- iso_files_areas()
    res <- cbind(Names, Areas)
    mzrange <- paste('mz:', round(iso_ref_eic()$mzs[input$iso_target_isotope, 1], 3), 
                     '-', round(iso_ref_eic()$mzs[input$iso_target_isotope, 2], 3))
    colnames(res) <- c('Names', mzrange) 
    res
  })
  
  quant_res <- reactiveVal(data.frame())
  
  quant_res_this <- reactive({
    if (input$iso_target_select == 'YES') {
      target <- input$iso_formula
      type <- input$iso_adduct
    } else {
      target <- input$iso_fmz
      type <- NA
    }
    rtmin <- input$iso_target_left
    rtmax <- input$iso_target_right
    mzmin <- round(iso_ref_eic()$mzs[input$iso_target_isotope, 1], 3)
    mzmax <- round(iso_ref_eic()$mzs[input$iso_target_isotope, 2], 3)
    areas <- iso_files_areas()
    this <- cbind(target, type, rtmin, rtmax, mzmin, mzmax, t(areas))
    colnames(this) <- c('target', 'type', 'rtmin', 'rtmax', 'mzmin', 'mzmax', iso_name_sample())
    
    this
  })
  
  observeEvent(input$add_button,{
    withProgress(message = 'Successful', value = 0.3, {
      quant_res(rbind(quant_res(), quant_res_this()))
    })
  })
  
  output$res_files_peaks <- renderTable({
    quant_res()
  })
  
  output$res_download <- downloadHandler(
    filename = "results.csv",
    content = function(file) {
      write.csv(quant_res(), file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
}