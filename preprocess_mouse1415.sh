#!/bin/bash

NETWORK_BASE="./networks/mouse1415/mouse1415"

#octave -q <<END
#s = load("1752-0509-4-140-s5.mat");
#model = s.mouse1415;
#[r2gi, r2gj, r2gv] = find(model.rxnGeneMat);
#
#a = [model.rxns(r2gi), model.genes(r2gj), num2cell(num2str(r2gv))];
#
#fid = fopen("./mouse1415.rxn2genes.tsv", "w");
#
#fprintf(fid, "rxn\tgene\tval\n");
#fprintf(fid, "%s\t%s\t%s\n", a'{:});
#
#fclose(fid);
#
#outputNetworkCytoscape(model, "mouse1415");
#END

#
#patch -p0 <<END
#--- mouse1415.sif       2013-08-14 12:15:44.706281184 +0400
#+++ mouse1415.sif       2013-08-14 12:15:58.442235920 +0400
#@@ -379,8 +379,8 @@
# EX_retinol-cis-11(e) rev retinol-cis-11(e)
# EX_retn(e) rev retn(e)
# EX_retnglc(e) rev retnglc(e)
#-EX_retpalm rev retpalm [deleted 10/09/2005  06:18:49 PM](e)
#-EX_retpalm(e) rev retpalm [deleted 10/09/2005  06:18:49 PM](e)
#+#EX_retpalm rev retpalm [deleted 10/09/2005  06:18:49 PM](e)
#+#EX_retpalm(e) rev retpalm [deleted 10/09/2005  06:18:49 PM](e)
# EX_rib-D(e) rev rib-D(e)
# EX_ribflv(e) rev ribflv(e)
# EX_s2l2fn2m2masn(e) rev s2l2fn2m2masn(e)
#--- mouse1415_nodeType.noa      2013-08-14 12:17:10.585998180 +0400
#+++ mouse1415_nodeType.noa      2013-08-14 12:17:15.937980540 +0400
#@@ -6466,7 +6466,7 @@
# retinol-cis-11(e) = met
# retn(e) = met
# retnglc(e) = met
#-retpalm [deleted 10/09/2005  06:18:49 PM](e) = met
#+#retpalm [deleted 10/09/2005  06:18:49 PM](e) = met
# rib-D(e) = met
# ribflv(e) = met
# s2l2fn2m2masn(e) = met
#--- mouse1415_nodeComp.noa	2013-08-29 19:08:39.062190515 -0500
#+++ mouse1415_nodeComp.noa	2013-08-29 19:08:49.062190970 -0500
#@@ -1335,7 +1335,7 @@
# retinol-cis-11(e) = e
# retn(e) = e
# retnglc(e) = e
#-retpalm [deleted 10/09/2005  06:18:49 PM](e) = deleted 10/09/2005  06:18:49 PM
#+#retpalm [deleted 10/09/2005  06:18:49 PM](e) = deleted 10/09/2005  06:18:49 PM
# rib-D(e) = e
# ribflv(e) = e
# s2l2fn2m2masn(e) = e
#@@ -3151,8 +3151,8 @@
# EX_retinol-cis-11(e) = e
# EX_retn(e) = e
# EX_retnglc(e) = e
#-EX_retpalm = deleted 10/09/2005  06:18:49 PM
#-EX_retpalm(e) = deleted 10/09/2005  06:18:49 PM
#+#EX_retpalm = deleted 10/09/2005  06:18:49 PM
#+#EX_retpalm(e) = deleted 10/09/2005  06:18:49 PM
# EX_rib-D(e) = e
# EX_ribflv(e) = e
# EX_s2l2fn2m2masn(e) = e
#END

#mv "$NETWORK_BASE_nodeComp."{noa,NA}
#mv "$WORKDIR/mouse1415_nodeType."{noa,NA}
#mv "$WORKDIR/mouse1415_edgeType."{noa,EA}
#rm "$WORKDIR/mouse1415_subSys.noa"

escape_primes() {
    f="$1"
    t=`mktemp`
    sed "s/'/_prime_/g" <"$f"  >"$t"
    mv "$t" "$f"
}

#escape_primes "$WORKDIR/mouse1415.rxn2genes.tsv"

./R/preprocess_network.R \
    -i $NETWORK_BASE \
    --mask "./misc/mets2mask.noncomp.hmdb.txt" \
    --cobra2hmdb "./misc/HMDB2cobraID_dictionary.noncomp.txt" \
    --metabolite-ids "./misc/metabolite_ids.tsv"
