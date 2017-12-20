# TarMet
TarMet is a shiny application for targeted metabolic analyses based on mass spectrometry.

## Installation  
<<<<<<< HEAD
### Install Depends: 
	install.packages(c("devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet"))
	library(devtools)
	httr::set_config( httr::config( ssl_verifypeer = 0L ) )
	install_github("zmzhang/baselineWavelet")
### Install TarMet
	install_github("hcji/TarMet")
	
=======
### R users:

	install.packages(c("devtools", "data.table", "enviPat", "Matrix", "shiny", "plotly"))
	source("https://bioconductor.org/biocLite.R")
    biocLite(c("MassSpecWavelet", "mzR"))
	library(devtools)
	httr::set_config( httr::config( ssl_verifypeer = 0L ) )
	install_github("zmzhang/baselineWavelet")

	install_github("hcji/TarMet")
	
### Not R users
download TarMet at [url](https://pan.baidu.com/s/1dEMfUF3) and install like a normal windows application.

>>>>>>> 35ddb0fd4d8068e97d21c71dc944db524b6ea1e7
## Usage:
  A user guide is included in the package	

## Contact
  For any questions, please contact:  ji.hongchao@foxmail.com