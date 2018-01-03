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
        fileInput('file', 'Choose one of your raw MS or Swtch-MS data files'),
        
        h4('Metabolite Information'),
        selectInput('target_select', 'How to define the targeted metabolite ?', c('formula', 'm/z of ion')),
        uiOutput('formula_contral'),
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
  ),
  
  tabPanel(
    'Swath-MS Analysis',
    titlePanel("Swath-MS Analysis"),
    sidebarLayout(
      sidebarPanel(
        numericInput('int_thres', 'Input the threshold value of peak intensity', 0.7, min=0),
        numericInput('corr_thres', 'Input the threshold value of coefficient of peak profile (from 0 to 1)', 0.7, min=0, max=1),
        numericInput('area_thres', 'Input the threshold value of coefficient of peak areas (from 0 to 1)', 0.7, min=0, max=1),
        numericInput('swath_ind', 'Input the sample index for plot', 1, step=1)
      ),
      mainPanel(
        h4('Fragment Ion EICs'),
        plotlyOutput('SwathEicPlot'),
        h4('Fragment Ion Peak Areas'),
        tableOutput('SwathAreaTable'),
        h4('Pseudo MS2'),
        tableOutput('SwathMS2')
      )
    )
  )
  
)