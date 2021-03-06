#!/usr/bin/env python
import re

with open("./glycan") as f:
    text = "".join(f.readlines())

descriptions = text.split("///\n")

m2name = ["met\tname\tformula\tpathway"]

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

    met_id = d["ENTRY"].split()[0]

    escaped_name = ""
    formula = ""
    if "NAME" in d:
        escaped_name = d["NAME"].replace("\n", " ") \
                                .replace('"', '\\"') \
                                .replace("<","") \
                                .replace(">","")
        if "FORMULA" in d:
            formula = d["FORMULA"]


    pathways = ""
    if "PATHWAY" in d:
        pathways = "+".join([p.split()[0] for p in d["PATHWAY"].split("\n")])

    m2name.append('%s\t"%s"\t"%s"\t"%s"' % (met_id, escaped_name, formula, pathways))

with open("gly2name.tsv", "w") as f:
    f.write("%s\n" % "\n".join(m2name))
