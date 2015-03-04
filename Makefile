GAM: description document examples build

GAM.db-package:
	make -C GAM.db

GAM.networks-package:
	make -C GAM.networks

r-dependencies:
	Rscript install_dependencies.R

fast-build: description
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
	mkdir -p GAM/data
	mv examplesGAM.rda GAM/data/

check: document
	R CMD check GAM
