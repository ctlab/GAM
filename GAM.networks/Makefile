build: document data/kegg.human.network.rda data/kegg.mouse.network.rda
	R -e 'library(devtools); build(".")'

fast-build:
	R -e 'library(devtools); build(".", vignettes=F)'

document:
	R -e 'library(devtools); document(".")'

data/kegg.human.network.rda data/kegg.mouse.network.rda data/kegg.yeast.network.rda data/kegg.arabidopsis.network.rda:
	mkdir -p data
	./make_networks.R
