GAM="../../GAM"

kegg.mouse.network.rda kegg.human.network.rda: \
		kegg.db.rda \
		make_networks.R
	./make_networks.R

kegg.db.rda: \
		enz2gene.tsv \
		rxn2enz.tsv \
		net.sif \
		rxn2name.tsv \
		met2name.tsv \
		mets2mask.lst \
		rxns2mask.lst \
		mets2collapse.tsv \
		make_kegg.db.R
	./make_kegg.db.R


met2name.tsv: cpd2name.tsv gly2name.tsv
	tail -n +2 gly2name.tsv | cat cpd2name.tsv - | tr -d "<>" > met2name.tsv

mets2collapse.tsv cpd2name.tsv: compound convert_compound.py
	./convert_compound.py

enz2gene.tsv: enzyme convert_enzyme.py
	./convert_enzyme.py

gly2name.tsv: glycan convert_glycan.py
	./convert_glycan.py

net.sif rxn2enz.tsv rxn2name.tsv: reaction convert_reaction.py
	./convert_reaction.py

rxns2mask.lst: reaction
	grep -B 6 "of \w*-step" reaction | grep ENTRY | grep -Po "R\d{5}" > rxns2mask.lst

ligand:
	wget -nH --cut-dirs=1 -r ftp://ftp.bioinformatics.jp/kegg/ligand

glycan: ligand ligand/glycan.tar.gz
	tar xzf ligand/glycan.tar.gz -C x glycan/glycan -O > glycan

compound: ligand ligand/compound.tar.gz
	tar xzf ligand/compound.tar.gz -C x compound/compound -O > compound

enzyme: ligand ligand/enzyme.tar.gz
	tar xzf ligand/enzyme.tar.gz -C x enzyme/enzyme -O > enzyme

reaction: ligand ligand/reaction.tar.gz
	tar xzf ligand/reaction.tar.gz -C x reaction/reaction -O > reaction

install: \
		kegg.mouse.network.rda \
		kegg.human.network.rda \
		kegg.db.rda
	cp kegg.*.rda ${GAM}/data/
