build: document data/kegg.db.rda data/met.id.map.rda
	R -e 'library(devtools); build(".")'

fast-build:
	R -e 'library(devtools); build(".", vignettes=F)'

document:
	R -e 'library(devtools); document(".")'

kegg/kegg.db.rda:
	make -C kegg kegg.db.rda

data/kegg.db.rda: kegg/kegg.db.rda
	mkdir -p data
	cp kegg/kegg.db.rda data/

metabocards/met.id.map.rda:
	make -C metabocards met.id.map.rda

data/met.id.map.rda: metabocards/met.id.map.rda
	mkdir -p data
	cp metabocards/met.id.map.rda data/
