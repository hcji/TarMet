# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Release note
v1.1.1 version is pre-released as the development version.

### What's new
1. New interface.
2. Allow load/save targeted compounds information via a config file.
3. Plot stack-bar figure automatically.
4. Fix a bug, which makes return error when no peaks detected.

### Known bug
1. Show error massage when re-input some parameters, but not affect the result.

### In the future
1. Allow run in command line.
2. Allow multi isotope tracer, like 13C and 15N.
3. DIA support.

## Installation  
### Install Depends
First, please install the depends with the following codes

	install.packages(c("rcdk", "tidyr", "devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet", "mzR"))

### Released version (Suggested)
Download the source package at [url](https://github.com/hcji/TarMet/releases/download/v1.1.0/TarMet_1.1.0.tar.gz) and install the package locally.

### Development version
Note: Development may include bugs/error and functions not ready.

	library(devtools)
	install_github("hcji/TarMet")

## Usage:

	library(TarMet)
	runTarMet()
	
  A user guide is included in the package, [Here](https://github.com/hcji/TarMet/releases/download/v1.1.1/TarMet.gif) is a gif of how to use the software.

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com
