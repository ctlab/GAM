all: description document build

fast-build: description
	R -e 'library(devtools); build("GAM", vignettes=F)'

build: description document
	R -e 'library(devtools); build("GAM")'

document:
	R -e 'library(devtools); document("GAM")'

description:
	./update_description.sh

examples:
	Rscript GAM/inst/make_examples.R
	mv examplesGAM.rda GAM/data/

check: document
	R CMD check GAM
