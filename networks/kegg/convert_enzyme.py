#!/usr/bin/env python
import re

reactions_file = "./enzyme"

with open(reactions_file) as f:
    text = "".join(f.readlines())

descriptions = text.split("///\n")

organism = "MMU"

e2g = ["enz\tgene"]

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
    genes = filter(lambda s: s.startswith(organism), genes)
    assert len(genes) <= 1
    if len(genes) == 0:
        continue
    genes = genes[0]
    genes = genes[len("%s: " % organism):].split(" ")
    genes_entrez = map(lambda s: re.sub("(.*)\\(.*\\)", "\\1", s), genes)
    e2g.extend(["%s\t%s" % (enz_id, gene) for gene in genes_entrez])


with open("enz2gene.tsv", "w") as f:
    f.write("%s\n" % "\n".join(e2g))
