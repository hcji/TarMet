# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Installation  

	install.packages(c("tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet", "mzR"))
	library(devtools)
	install_github("zmzhang/baselineWavelet")

	install_github("hcji/TarMet")
	library(TarMet)
	runTarMet()

## Usage:
  A [user guide](https://github.com/hcji/TarMet/blob/master/inst/doc/TarMet.pdf) is included in the package	

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
