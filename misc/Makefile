met.id.map.rda: metabocards-2.5.extr.tsv ../networks/kegg/met2name.tsv preprocess_id_maps.R
	./preprocess_id_maps.R

metabocards-2.5.extr.tsv: metabocards-2.5.tmp.txt
	paste - - - < metabocards-2.5.tmp.txt | \
		sed "s/^#BEGIN_METABOCARD //" | \
		sed "s/Not available/NA/i" | \
		sed "s/^\(.*\)\t\(\w\+\); /\1\t\2\n\1\t/" | \
		sed "s/\(HMDB00656\t\w*\t\)NA/\1C05526/" | \
		sed "1i HMDB\tformula\tKEGG" > metabocards-2.5.extr.tsv

metabocards-2.5.tmp.txt: metabocards-2.5.zip
	unzip -p metabocards-2.5.zip | \
		grep -A 1 -e "^# kegg_compound_id:" -e "BEGIN_METABOCARD" -e "^# chemical_formula:" | \
		grep -v -e -- -e '^$$' -e "chemical_formula" -e "kegg_compound_id" > metabocards-2.5.tmp.txt


metabocards-2.5.zip:
	wget "http://www.hmdb.ca/downloads/2_5/metabocards.zip" -O metabocards-2.5.zip

#metabocards-1.0.txt.gz:
#	wget "http://www.hmdb.ca/downloads/1_0/metabocards.txt.gz" -O metabocards-1.0.txt.gz

install: met.id.map.rda
	cp met.id.map.rda ../GAM/data/

