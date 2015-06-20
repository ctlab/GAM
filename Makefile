GAM: description document examples build

GAM.db-package: r-dependencies
	make -C GAM.db

GAM.networks-package: r-dependencies
	make -C GAM.networks

r-dependencies:
	Rscript install-dependencies.R

fast-build: description r-dependencies
	R -e 'library(devtools); build("GAM", vignettes=F)'

build: description document r-dependencies
	R -e 'library(devtools); build("GAM")'

build.db: description document r-dependencies
	R -e 'library(devtools); build("GAM.db")'

document: r-dependencies
	R -e 'library(devtools); document("GAM")'

description:
	./update_description.sh

examples: r-dependencies
	Rscript GAM/inst/make_examples.R
	mkdir -p GAM/data
	mv examplesGAM.rda GAM/data/

check: document
	R CMD check GAM

test:
	R -e 'library("devtools"); library("testthat"); load_all("GAM"); test_dir("GAM/inst/tests/")'
