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

mets2collapse.tsv cpd2name.tsv: compound process_compound.py
	./process_compound.py

enz2gene.tsv: enzyme process_enzyme.py
	./process_enzyme.py

gly2name.tsv: glycan process_glycan.py
	./process_glycan.py

net.sif rxn2enz.tsv rxn2name.tsv: reaction process_reaction.py
	./process_reaction.py

rxns2mask.lst: reaction
	grep -B 6 "of \w*-step" reaction | grep ENTRY | grep -Po "R\d{5}" > rxns2mask.lst

ligand:
	wget -nH --cut-dirs=1 -r ftp://ftp.bioinformatics.jp/kegg/ligand

glycan: ligand ligand/glycan.tar.gz
	tar xzf ligand/glycan.tar.gz glycan/glycan -O > glycan

compound: ligand ligand/compound.tar.gz
	tar xzf ligand/compound.tar.gz compound/compound -O > compound

enzyme: ligand ligand/enzyme.tar.gz
	tar xzf ligand/enzyme.tar.gz enzyme/enzyme -O > enzyme

reaction: ligand ligand/reaction.tar.gz
	tar xzf ligand/reaction.tar.gz reaction/reaction -O > reaction
