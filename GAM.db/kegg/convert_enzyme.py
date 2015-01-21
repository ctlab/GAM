#!/usr/bin/env python
import re

with open("./enzyme") as f:
    text = "".join(f.readlines())

descriptions = text.split("///\n")

organisms = ["MMU", "HSA"]

e2g = ["enz\tgene\torganism"]
e2g_all = ["enz\tgene\torganism"]

for description in descriptions:
    if description == "":
        continue
    lines = description.split("\n")
    shift = 12
    last_name = None
    d = {}
    for line in lines:
        if line == "":
            continue
        name = line[:shift].strip()
        value = line[shift:]
        if name == "":
            name = last_name

        if name in d:
            d[name] += "\n" + value
        else:
            d[name] = value

        last_name = name

    enz_id = d["ENTRY"].split()[1]

    if not "GENES" in d:
        d["GENES"] = ""

    genes = d["GENES"].split("\n")
    for organism_genes in genes:
        if len(organism_genes) == 0:
            continue

        (organism, organism_genes) = organism_genes.split(":", 1)
        organism_genes = organism_genes.split(" ")
        genes_entrez = map(lambda s: re.sub("(.*)\\(.*\\)", "\\1", s), organism_genes)
        if organism in organisms:
            e2g.extend(["%s\t%s\t%s" % (enz_id, gene, organism) for gene in genes_entrez if gene != ""])
        e2g_all.extend(["%s\t%s\t%s" % (enz_id, gene, organism) for gene in genes_entrez if gene != ""])


with open("enz2gene.tsv", "w") as f:
    f.write("%s\n" % "\n".join(e2g))

with open("enz2gene.all.tsv", "w") as f:
    f.write("%s\n" % "\n".join(e2g_all))

