#!/usr/bin/env python
import re

reactions_file = "./compound"

with open(reactions_file) as f:
    text = "".join(f.readlines())

descriptions = text.split("///\n")

m2name = ["met\tname"]

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

    if "NAME" in d:
        m2name.append('%s\t"%s"' % (met_id, d["NAME"].replace("\n", " ")))


with open("met2name.tsv", "w") as f:
    f.write("%s\n" % "\n".join(m2name))
