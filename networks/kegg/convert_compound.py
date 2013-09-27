#!/usr/bin/env python
import re

with open("./compound") as f:
    text = "".join(f.readlines())

descriptions = text.split("///\n")

m2name = ["met\tname\tpathway"]

anomers = {}

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
    if "NAME" in d:
        escaped_name = d["NAME"].replace("\n", " ") \
                                .replace('"', '\\"') \
                                .replace("<","") \
                                .replace(">","")

        if "FORMULA" in d:
            for name in d["NAME"].split("\n"):
                if name.endswith(";"):
                    name = name[:-1]

                name_full = name

                if name == "beta-Tyrosine":
                    pass
                elif name.startswith("beta-"):
                    name = name[len("beta-"):]
                elif name.startswith("alpha-"):
                    name = name[len("alpha-"):]




                name = "%s: %s" % (d["FORMULA"], name)

                if not name in anomers:
                    anomers[name] = []
                anomers[name].append((name_full, met_id))


    pathways = ""
    if "PATHWAY" in d:
        pathways = "+".join([p.split()[0] for p in d["PATHWAY"].split("\n")])

    m2name.append('%s\t"%s"\t"%s"' % (met_id, escaped_name, pathways))


with open("mets2collapse.tsv", "w") as f:
    f.write('short\tfull\tfrom\tto\n')
    for (name, mets) in anomers.iteritems():
        if (len(set(zip(*mets)[1])) == 1):
            continue

        base_id = sorted(set(zip(*mets)[1]))[0]

        for (name_full, met_id) in mets:
            f.write('"%s"\t"%s"\t%s\t%s\n' % (name, name_full, met_id, base_id))


with open("cpd2name.tsv", "w") as f:
    f.write("%s\n" % "\n".join(m2name))
