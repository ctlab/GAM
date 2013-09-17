#!/usr/bin/env python
import re

reactions_file = "./compound"

with open(reactions_file) as f:
    text = "".join(f.readlines())

descriptions = text.split("///\n")

m2name = ["met\tname"]

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


    if "NAME" in d:
        m2name.append('%s\t"%s"' % (met_id, d["NAME"].replace("\n", " ")))
        if "FORMULA" in d:
            for name in d["NAME"].split("\n"):
                if name.endswith(";"):
                    name = name[:-1]

                name_full = name

                if name.startswith("beta-"):
                    name = name[len("beta-"):]
                elif name.startswith("alpha-"):
                    name = name[len("alpha-"):]




                name = "%s: %s" % (d["FORMULA"], name)

                if d["FORMULA"] == "C6H13O9P":
                    print name, met_id, name_full

                if not name in anomers:
                    anomers[name] = []
                anomers[name].append((name_full, met_id))



with open("mets2collapse.tsv", "w") as f:
    f.write('short\tfull\tfrom\tto\n')
    for (name, mets) in anomers.iteritems():
        if (len(set(zip(*mets)[1])) == 1):
            continue

        base_id = sorted(set(zip(*mets)[1]))[0]

        for (name_full, met_id) in mets:
            f.write('"%s"\t"%s"\t%s\t%s\n' % (name, name_full, met_id, base_id))


#with open("met2name.tsv", "w") as f:
#    f.write("%s\n" % "\n".join(m2name))
