# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Installation  
### Install Depends: 
	install.packages(c("devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet"))
	library(devtools)
	httr::set_config( httr::config( ssl_verifypeer = 0L ) )
	install_github("zmzhang/baselineWavelet")
### Install TarMet
	install_github("hcji/TarMet")
	
## Usage:
  A user guide is included in the package	

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com