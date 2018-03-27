library(shiny)
library(plotly)

data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

shinyUI(fluidPage(
  titlePanel(img(src="logo.png", width="20%"), "TarMet"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 fileInput('files', 'Choose multiple raw data files', multiple=TRUE),
                 numericInput('resolution', 'Input the resolution of your MS instrument', 50000),
                 uiOutput('sampleCtrl'),
                 selectInput('input', 'Select how to input the target compounds', c('config file', 'direct')),
                 selectInput('type', 'Select the type of the assay', c('targeted analysis', 'isotopic tracer')),
                 selectInput('define', 'Select how to define the target compound', c('formula','mass-to-charge')),
                 selectInput('adduct', 'Select the adduct type', choices = list(
                   Positive = adducts$Name[adducts$Ion_mode == 'positive'],
                   Negative = adducts$Name[adducts$Ion_mode == 'negative']
                 )),
                 uiOutput('targetCtrl1'),
                 uiOutput('targetCtrl2'),
                 uiOutput('formulaCtrl'),
                 uiOutput('tracerCtrl1'),
                 uiOutput('tracerCtrl2'),
                 actionButton('confirm', 'Confirm')
                 ),
                 
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel('Isotopologues Analysis',
                           sidebarLayout(
                             sidebarPanel(width = 3,
                                          h4('Peak Detection'),
                                          selectInput('baseline', 'Remove Baseline?', c(TRUE, FALSE)),
                                          selectInput('smooth', 'Smooth EIC?', c(FALSE, TRUE)),
                                          numericInput('snr.th', 'SNR threshold', 5),
                                          numericInput('scale.th', 'Scale threshold (s)', 5),
                                          numericInput('int.th', 'Intensity threshold (above the baseline)', 0),
                                          uiOutput('targetRtCtrl')
                                          ),
                             mainPanel(
                               plotlyOutput('EICPlot'),
                               tableOutput('PeakInfo'),
                               tableOutput('PeakArea')
                               
                             )
                           )),
                  tabPanel('Quantitative Analysis'),
                  tabPanel('Result List')
    )
    
  ))
))
