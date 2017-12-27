# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Installation  

### R users:

	install.packages(c("tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet", "mzR"))
	library(devtools)
	install_github("zmzhang/baselineWavelet")

	install_github("hcji/TarMet")
	library(TarMet)
	runTarMet()
	
### Not R users
Download TarMet at [url](https://www.researchgate.net/publication/322065923_Setup_file_of_TarMet_0991_version) and install like a normal windows application.
Note, it may take several minutes to install packages when starts for the first time.

## Usage:
  A [user guide](https://github.com/hcji/TarMet/blob/master/inst/doc/TarMet.pdf) is included in the package	

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
