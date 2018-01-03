# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Installation  

### Released version (suggested)

Download the source package at [url](https://github.com/hcji/TarMet/archive/1.0.1.tar.gz) and install the package locally.

### Development version

	install.packages(c("tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet", "mzR"))
	install_github("hcji/TarMet")

## Usage:

	library(TarMet)
	runTarMet()
	
  A [user guide](https://github.com/hcji/TarMet/blob/master/vignettes/TarMet.Rmd) is included in the package	

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
