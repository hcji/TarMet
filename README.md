## TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.


### In the future
1. Support data-independent acquisition.
2. Support MS/MS library verification.

### Installation  

### Released version

Download the source package at [url](https://github.com/hcji/TarMet/releases/download/v1.1.1/TarMet_1.1.1.tar.gz) and install the package locally.

### Development version

	install.packages(c("BiocManager", "rcdk", "tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	BiocManager::install("MassSpecWavelet")
	BiocManager::install("mzR")
	
	library(devtools)
	install_github("hcji/TarMet")

### Usage:

	library(TarMet)
	runTarMet()
	
  A user guide is included in the package, [Here](https://github.com/hcji/TarMet/releases/download/v1.1.1/TarMet.gif) is a gif of how to use the software.

### Command line:

  Detailed description of using TarMet with R scripts can be found [Here](https://github.com/hcji/TarMet/blob/master/vignettes/TarMetCL.Rmd) and [Here](https://rpubs.com/jihongchao/TarMetCL).
  
### Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
