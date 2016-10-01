#!/usr/bin/env python
import re
import sys
from itertools import *
from pprint import pprint
from collections import deque, Counter

import logging
from logging import debug, info, warn, error

logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)-15s - %(levelname)s: %(message)s"
        )

align_re = re.compile("\\d+ +(\\d+):\\w+ +(\\d+):\\w+")
def get_align_numbers(s):
    m = align_re.match(s)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    else:
        return (None, None)


def append_entry(d, name, value):
    if name in d:
        if isinstance(d[name], list):
            d[name].append(value)
        elif not d[name].endswith(" "):
            d[name] = [d[name], value]
        else:
            d[name] += value
    else:
        d[name] = value

def iter_descriptions(lines):
    shift = 12

    d = {}
    last_name = None


    for line in lines:
        line = line.rstrip("\n")
        if line == "":
            continue

        if line == "///":
            d["ENTRY"] = d["ENTRY"].split()[0]
            yield d
            d = {}

        name = line[:shift].strip()
        value = line[shift:]

        if line.startswith(" "):
            if not supentry and name != "":
                supentry = last_name
                d[supentry] = {}
        else:
            supentry = None

        if name == "":
            name = last_name

        if supentry is None:
            append_entry(d, name, value)
        else:
            append_entry(d[supentry], name, value)

        last_name = name



class Atom:
    def __init__(self, i, element, type, bonds=None, comment=""):
        self.i = i
        self.element = element
        self.type = type
        if bonds is None:
            bonds = set()
        self.bonds = bonds
        pass

    def __str__(self):
        return "%d:%s" % (self.i, self.type)

def parse_atom(s):
    ss = s.split()
    (i, type, element, x, y) = ss[:5]
    comment = ss[5] if len(ss) > 5 else ""
    return Atom(int(i) - 1, element, type, comment=comment)

class Bond:
    def __init__(self, i, start, end, valence, comment=""):
        self.i = i
        self.start = start
        self.end = end
        self.valence = valence
        self.comment = comment

    def opposite(self, start):
        return self.end if self.start == start else self.start

    def __str__(self):
        return "%d:%d-%s-%d" % (self.i, self.start, self.valence, self.end)

def parse_bond(s):
    ss = s.split()
    (i, start, end, valence) = ss[:4]
    comment = ss[4] if len(ss) > 4 else ""
    return Bond(int(i) - 1, int(start) - 1, int(end) - 1, valence, comment)

class Structure:
    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds

        for (i, bond) in izip(count(), self.bonds):
            self.atoms[bond.start].bonds.add(bond)
            self.atoms[bond.end].bonds.add(bond)

    def _plain(self):
        return (
                [ a.element for a in self.atoms ],
                [ (b.start, b.end) for b in self.bonds ]
                )

    def __eq__(self, other):
        if isinstance(other, Structure):
            return self._plain() == other._plain()
        return NotImplemented

    def __str__(self):
        return "Structure(%d, %d)" % (len(self.atoms), len(self.bonds))

def parse_structure(atoms_s, bonds_s):
    atoms = map(parse_atom, atoms_s[1:])
    bonds = map(parse_bond, bonds_s[1:])
    return Structure(atoms, bonds)

compound_structure = {}

with open("compound", "r") as compound_in:
    for d in iter_descriptions(compound_in):
        compound_id = d['ENTRY']
        if "ATOM" in d and "BOND" in d:
            compound_structure[compound_id] = parse_structure(d['ATOM'], d['BOND'])


def remaps(st1, st2, filter=['C', 'N']):
    n = len(st1.atoms)
    m = len(st1.bonds)
    assert len(st2.atoms) == n
    assert len(st2.bonds) == m

    good_maps = []

    if n == 0:
        return good_maps

    marked = [False] * n


    cur_queue = range(n)
    #def dfs(u):
    #    if u in cur_queue:
    #        return

    #    cur_queue.append(u)

    #    at1 = st1.atoms[u]

    #    for b1 in at1.bonds:

    #        v = b1.opposite(u)

    #        dfs(v)

    #dfs(0)

    assert len(cur_queue) == n

    cur_map = [None] * n
    cur_map_reverse = [None] * n

    filtered_atoms = [ i for i in xrange(n) if st1.atoms[i].element in filter]
    filtered_good_maps = []

    if st1 == st2:
        return dict([(i, i) for i in filtered_atoms])

    common_map = []

    def remaps_impl(i = 0):
        if len(filtered_good_maps) >= 4:
            common_map[0] = {}
            return

        if len(common_map) == 1 and len(common_map[0]) == 0:
            return

        if i == n:
            filtered_gm = dict([(i, cur_map[i]) for i in filtered_atoms])
            if not filtered_gm in filtered_good_maps:
                filtered_good_maps.append(filtered_gm)

            if len(common_map) == 0:
                common_map.append(filtered_gm)

            common_map[0] = dict(common_map[0].viewitems() & filtered_gm.viewitems())

            return

        u = cur_queue[i]
        at1 = st1.atoms[u]

        candidates = set([j for j in xrange(n) if (st2.atoms[j].type == at1.type)])
        #print "precandidates for %s:" % at1
        #print([str(st2.atoms[j]) for j in candidates])
        candidates = set([j for j in xrange(n) if (st2.atoms[j].element == at1.element) and (cur_map_reverse[j] is None)])
        #print "candidates for %s:" % at1
        #print([str(st2.atoms[j]) for j in candidates])

        for b1 in at1.bonds:
            #print "    ", str(b1)
            v = b1.opposite(u)
            
            if cur_map[v] is None:
                continue

            v2 = cur_map[v]
            atx2 = st2.atoms[v2]
            #print "    ", str(atx2)
            #print "    ", [str(b2) for b2 in atx2.bonds]


            candidates1 = set([b2.opposite(v2) for b2 in atx2.bonds if \
                    True])
            #print "    ", [str(st2.atoms[j]) for j in candidates1]
            candidates1 = set([b2.opposite(v2) for b2 in atx2.bonds if \
            #        b2.valence == b1.valence and \
                    st2.atoms[b2.opposite(v2)].element == at1.element])
            #print "    ", [str(st2.atoms[j]) for j in candidates1]

            candidates = candidates.intersection(candidates1)
            #print [str(st2.atoms[j]) for j in candidates]

        #print "#########"


        for j in candidates:
            #print "mapping  %s -> %s" % (str(at1), str(st2.atoms[j]))
            cur_map[u] = j
            cur_map_reverse[j] = u
            remaps_impl(i + 1)
            cur_map_reverse[j] = None
            #print "---- unmapping  %s -> %s" % (str(at1), str(st2.atoms[j]))


        cur_map[u] = None
        return

    remaps_impl()

    #for good_map in filtered_good_maps:
    #    pprint(["%s <-> %s" % (st1.atoms[i], st2.atoms[good_map[i]]) for i in filtered_atoms])
    if len(common_map) == 0:
        return None
    return common_map[0]


def structure_ok(entry):
    met = entry['COMPOUND']
    return met in compound_atom and compound_atom[met] == strip_coordinates(entry['ATOM']) and \
           met in compound_bond and compound_bond[met] == entry['BOND']

bad_compounds = set([
        # phytate and its dirivatives
        "C01204", "C04563", "C04579", "C01284", "C11174", "C11526", "C15990", "C15991"
        ])

def normalize_entry(d):
    met = d['COMPOUND']
    if met in bad_compounds:
        warn("Skipping %s because bad" % met)
        return None

    entry_st = parse_structure(d['ATOM'], d['BOND'])

    if not met in compound_structure:
        return None

    ref_st = compound_structure[met]

    if len(ref_st.atoms) >= 100:
        warn("Skipping %s because >= 100 atoms" % met)
        return None

    if len(ref_st.atoms) != len(entry_st.atoms):
        warn("Different number of atoms for %s" % met)
        return None

    if len(ref_st.bonds) != len(entry_st.bonds):
        warn("Different number of bonds for %s" % met)
        return None

    if Counter([a.element for a in ref_st.atoms]) != \
       Counter([a.element for a in entry_st.atoms]):
        warn("Different sets of atoms for %s" % met)
        return None

    m = remaps(entry_st, ref_st)

    if m is None:
        warn("No remaps for %s" % met)
        return None
    assert isinstance(m, dict), m

    expect_maps = len([a for a in entry_st.atoms if a.element in ['C', 'N']])

    if len(m) == 0 and expect_maps > 0:
        warn("Too many remaps for %s" % met)
        return None


    if len(m) < expect_maps:
        warn("Not all atoms mapped for %s" % met)

    return m


rpaligns = ["\t".join(["rpair", "atom.x", "atom.y"])]

print "processing RPAIRs"

skipped = []
with open("rpair", "r") as rpair_in:
    for d in iter_descriptions(rpair_in):

        rpair_id = d['ENTRY']

        #if rpair_id != "RP02527":
        #    continue

        print rpair_id
        met1 = d['ENTRY1']['COMPOUND']
        met2 = d['ENTRY2']['COMPOUND']

        met1_map = normalize_entry(d['ENTRY1'])
        if met1_map is None:
            skipped.append((rpair_id, met1, met2))
            continue

        met2_map = normalize_entry(d['ENTRY2'])

        if met2_map is None:
            skipped.append((rpair_id, met1, met2))
            continue

        aligns = d["ALIGN"][1:]

        # there can be maps C -> R
        rpaligns.extend(
                ["%s\t%s_%s\t%s_%s" % (rpair_id, met1, met1_map[a1 - 1],
                                                 met2, met2_map[a2 - 1]) 
                    for (a1, a2) in map(get_align_numbers, aligns)
                    if (not a1 is None) and ((a1-1) in met1_map) \
                                        and ((a2-1) in met2_map)])

with open("skipped.tsv", "w") as f:
    for s in skipped:
        f.write("%s\t%s\t%s\n" % s)


with open("rpaligns.tsv", "w") as f:
    f.write("%s\n" % "\n".join(rpaligns))
