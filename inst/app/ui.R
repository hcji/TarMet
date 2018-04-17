library(shiny)
library(plotly)

data("isotopes", package = "enviPat")
data("adducts", package = "enviPat")

shinyUI(fluidPage(
  titlePanel(img(src="logo.png", width="20%"), "TarMet"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 fileInput('files', 'Choose multiple raw data files', multiple=TRUE),
                 numericInput('resolution', 'Resolution of your MS instrument:', 50000),
                 uiOutput('sampleCtrl'),
                 selectInput('input', 'How to input the target compounds?', c('config file', 'new compound')),
                 selectInput('type', 'Type of the assay:', c('targeted analysis', 'isotopic tracer', 'data independent analysis'))
                 ),
                 
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel('Isotopologues Analysis',
                           sidebarLayout(
                             sidebarPanel(width = 3,
                                          h4('Definition'),
                                          uiOutput('targetCtrl1'),
                                          uiOutput('targetCtrl2'),
                                          uiOutput('diaCtrl'),
                                          uiOutput('tracerCtrl1'),
                                          uiOutput('tracerCtrl2'),
                                          uiOutput('paraCtrl'),
                                          actionButton('confirm', 'Confirm'),
                                          actionButton('addConfig', 'Add to Config'),
                                          uiOutput('targetRtCtrl')
                                          
                                          ),
                             mainPanel(
                               plotlyOutput('EICPlot'),
                               tableOutput('PeakInfo'),
                               tableOutput('PeakArea')
                             )
                           )),
                  
                  tabPanel('Quantitative Analysis',
                           sidebarLayout(
                             sidebarPanel(width = 3,
                                          selectInput('alignment', 'Alignment EICs across samples?', c(TRUE, FALSE)),
                                          uiOutput('alignmentCtrl'),
                                          uiOutput('isoPlotCtrl'),
                                          actionButton('addButton', 'Add to List')
                                          ),
                             mainPanel(
                               plotlyOutput('SamplesPlot'),
                               tableOutput('SamplesArea'),
                               plotlyOutput('SamplesBarPlot')
                             )
                           )),
                             
                  tabPanel('Result List',
                           h3('Result List'),
                           tableOutput('resultList'),
                           downloadButton("resultDown", "Download"),
                           h3('Config List'),
                           tableOutput('configList'),
                           downloadButton("configDown", "Download")
                           )
    )
    
  ))
))
