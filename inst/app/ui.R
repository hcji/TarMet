data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

navbarPage(
  "Targeted Bio-Mass Analyzer",
  
  tags$head(
    tags$link(
      rel = "icon", 
      type = "image/x-icon", 
      href = "http://localhost:1984/default.ico")
  ),
  
  tabPanel(
    'Isotopic Analysis',
    titlePanel("Isotopic Analysis"),
    sidebarLayout(
      sidebarPanel(
        fileInput('file', 'Choose one of your raw data files'),
        
        h4('Metabolite Information'),
        textInput('formula', 'Input the targeted metabolite'),
        numericInput('fmz', 'Input the monoisotopic mass of the targeted metabolite (if the formula is unknow)', -1),
        selectInput('adduct', 'Select the type of adduct', choices = list(
          Positive = adducts$Name[adducts$Ion_mode == 'positive'],
          Negative = adducts$Name[adducts$Ion_mode == 'negative']
        )),
        actionButton("button_formula", "Confirm"),
        
        h4('Isotopic Information (only used when the formula is given)'),
        numericInput('threshold', 'Input the threshold of the relative abundance of isotopic peaks', 0.001),
        numericInput('nmax', 'Input n, where at most M+n isotopologues are detected.', 3),
        
        h4('EIC Information'),
        numericInput('ppm', 'Input the tolerance of the m/z difference (ppm)', 100),
        uiOutput('rtControl'),
        selectInput('ifsmooth', 'Select wether the baselines are removed or not', c(TRUE, FALSE)),
        
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
        numericInput('align.shift', 'Input the maximum retention time of shift (s)', 20),
        numericInput('align.seg', 'Input the size of segment of PAFFT (s)', 20),
        
        h4('Isotope Information'),
        uiOutput('EICs_Control'),
        
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