# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Installation  

### Released version

Download the source package at [url](https://github.com/hcji/TarMet/releases) and install the package locally.

### Development version

	install.packages(c("rcdk", "tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet", "mzR"))
	
	library(devtools)
	install_github("hcji/TarMet")

## Usage:

	library(TarMet)
	runTarMet()
	
  A [user guide](https://github.com/hcji/TarMet/blob/master/vignettes/TarMet.Rmd) is included in the package	

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
