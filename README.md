# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Note
The source code roll back to v1.1.1 version, which is identical to the latest release version. The support of DIA dataset is removed. 
New package for DIA-MS is developping. Pay close attention to [TarSWATH](https://github.com/hcji/TarSWATH).

### Known bug
1. Show error massage when re-input some parameters, but not affect the result.

### In the future
1. Allow run in command line.
2. Allow multi isotope tracer, like 13C and 15N.

## Installation  

### Released version

Download the source package at [url](https://github.com/hcji/TarMet/releases/download/v1.1.1/TarMet_1.1.1.tar.gz) and install the package locally.

### Development version

	install.packages(c("BiocManager", "rcdk", "tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	BiocManager::install("MassSpecWavelet")
	BiocManager::install("mzR")
	
	library(devtools)
	install_github("hcji/TarMet")

## Usage:

	library(TarMet)
	runTarMet()
	
  A user guide is included in the package, [Here](https://github.com/hcji/TarMet/releases/download/v1.1.1/TarMet.gif) is a gif of how to use the software.

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
