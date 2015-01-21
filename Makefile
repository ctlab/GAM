GAM: build

GAM.db:
	make -C GAM.db

GAM.networks:
	make -C GAM.networks

all: description document examples build

r-dependencies:
	Rscript install_dependencies.R

fast-build: description data
	R -e 'library(devtools); build("GAM", vignettes=F)'

build: description document
	R -e 'library(devtools); build("GAM")'

build.db: description document
	R -e 'library(devtools); build("GAM.db")'

document:
	R -e 'library(devtools); document("GAM")'

description:
	./update_description.sh

examples:
	Rscript GAM/inst/make_examples.R
	mv examplesGAM.rda GAM/data/

check: document
	R CMD check GAM
