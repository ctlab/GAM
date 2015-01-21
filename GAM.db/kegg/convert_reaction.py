#!/usr/bin/env python
import re

with open("./reaction") as f:
    text = "".join(f.readlines())

def get_compound(s):
    return s.split()[-1]

def compounds(s):
    s = re.sub("\(.*?\)", "", s)
    return map(get_compound, s.split("+"))

descriptions = text.split("///\n")

m2m = ["met.x\trxn\tmet.y"]
r2e = ["rxn\tenz"]

r2name = ["rxn\tname\tpathway"]
rpairs = ["\t".join(["rxn", "rpair", "met.x", "met.y", "rptype"])]

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

    rxn_id = d["ENTRY"].split()[0]

    if not "ENZYME" in d:
        d["ENZYME"] = ""

    enzymes = d["ENZYME"].split()
    r2e.extend(["%s\t%s" % (rxn_id, enzyme) for enzyme in enzymes])

    equation = d["EQUATION"]

    (left, right) = map(compounds, equation.split("<=>"))
    m2m.extend(["%s\t%s\t%s" % (c1, rxn_id, c2) for c1 in left for c2 in right])
    escaped_name = ""
    if "NAME" in d:
        escaped_name = d["NAME"].replace("\n", " ") \
                                .replace('"', '\\"') \
                                .replace("<","") \
                                .replace(">","")


    pathways = ""
    if "PATHWAY" in d:
        pathways = "+".join([p.split()[0] for p in d["PATHWAY"].split("\n")])

    r2name.append('%s\t"%s"\t"%s"' % (rxn_id, escaped_name, pathways))

    if "RPAIR" in d:
        for rpair in d["RPAIR"].split("\n"):
            if len(rpair) == 0:
                continue
            # :ToDo: grab reaction class too (4th field)
            (rpair, mets, rptype) = rpair.split()[0:3]
            (metx, mety) = mets.split("_")
            rpairs.append("%s\t%s\t%s\t%s\t%s" % (rxn_id, rpair, metx, mety, rptype))

with open("net.sif", "w") as f:
    f.write("%s\n" % "\n".join(m2m))

with open("rxn2enz.tsv", "w") as f:
    f.write("%s\n" % "\n".join(r2e))

with open("rxn2name.tsv", "w") as f:
    f.write("%s\n" % "\n".join(r2name))

with open("rpairs.tsv", "w") as f:
    f.write("%s\n" % "\n".join(rpairs))
