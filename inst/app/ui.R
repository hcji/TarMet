navbarPage(
  "Targeted Bio-Mass Analyzer",
  
  tabPanel(
    'Isotopologues Analysis',
    titlePanel("Isotopologues Analysis"),
    sidebarLayout(
      sidebarPanel(
        fileInput('iso_files', 'Choose multiple raw data files', multiple=TRUE),
        uiOutput('iso_ctrl_sample'),
        
        h4('Metabolite Information'),
        selectInput('iso_target_select', 'How to define the targeted metabolite ?', c('formula', 'm/z of ion')),
        uiOutput('iso_ctrl_target'),
        actionButton('iso_target_button', 'Confirm'),

        h4('EIC Information'),
        numericInput('iso_ppm', 'Input the tolerance of the m/z difference (ppm)', 50),
        uiOutput('iso_ctrl_rt'),
        selectInput('iso_baseline', 'Select wether the baselines are removed or not', c(TRUE, FALSE)),

        h4('Peak Detection Information'),
        numericInput('iso_peak_snr.th', 'Input the minimum SNR of peaks', 5),
        numericInput('iso_peak_scale.th', 'Input the scale range of the peak (s)', 5),
        numericInput('iso_peak_int.th', 'Input the minimal absolute intensity (above the baseline) of peaks to be picked', 0),

        h4('Target Peak Defination'),
        uiOutput('iso_ctrl_targetRT'),
        
        h4('Alignment Information'),
        selectInput('iso_align', 'Select wether the EICs are aligned across samples or not', c(TRUE, FALSE)),
        numericInput('iso_align.shift', 'Input the maximum retention time of shift (s)', 20),
        numericInput('iso_align.seg', 'Input the size of segment of PAFFT (s)', 20),
        
        h4('Confirm the Result'),
        actionButton('add_button', 'Add to List'),
        downloadButton("iso_download", "Download")
      ),
      
      mainPanel(
        h3('Isotopologues Analysis Result of Reference Sample'),
        h4('Isotopologues EICs'),
        plotlyOutput('iso_eic_plot'),
        h4('Peak Information'),
        tableOutput('iso_peak_info'),
        tableOutput('iso_peak_area'),
        h3('Quantitative Analysis Result of All Samples'),
        h4('EICs of targeted isotopologues of samples'),
        plotlyOutput('iso_files_EIC_Plot'),
        h4('Peak areas of targeted isotopologues of samples'),
        tableOutput('iso_files_peaks')
      )
    )
  ),
  
  tabPanel(
    'Quantitative Result',
    titlePanel("Quantitative Result"),
    tableOutput('res_files_peaks'),
    downloadButton("res_download", "Download")
  )
)