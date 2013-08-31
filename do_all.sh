#!/bin/sh

# setting exit on error
set -e


GENE_EXPRS="$1"
GENE_CONDITIONS="$2"
MET_EXPRS="$3"
MET_CONDITIONS="$4"
WORKDIR_PREFIX="$5"

export ILOG_LICENSE_FILE="$HOME/Dropbox/4alexey_fedor/Network_Analysis_CPLEX/access.ilm"


#NETWORK_BASE="networks/test/test"
NETWORK_BASE="networks/mouse1415/mouse1415"


echo "LPS_IFNG vs IL4 "
STATE1=MandLPSandIFNg
STATE2=MandIL4
FDR="1e-7"

#echo "IL4 vs MO"
#METS_BASE=Ctrl-MandIL4
#GENES_BASE=MandIL4-MO
#FDR="5e-3"

#echo "LPS_IFNG vs MO"
#METS_BASE=Ctrl-MandLPSandIFNg
#GENES_BASE=MO-MOandLPSandIFNg
#FDR="1e-4"

WORKDIR="${WORKDIR_PREFIX}_$STATE1-$STATE2"
mkdir $WORKDIR || true


MET_PVALS=$WORKDIR/Mets.pval
GENE_PVALS=$WORKDIR/Genes.pval
RXN_PVALS=$WORKDIR/Rxns.pval

MET_COBRA_PVALS=$WORKDIR/Mets.cobra.pval
MET_NONCOMP_PVALS=$WORKDIR/Mets.cobra.noncomp.pval

echo $STATE1 vs $STATE2

#echo "Preprocessing mouse1415 network"
#./preprocess_mouse1415.sh

#echo "Metabolite differential expression"
#./R/LogDiffExpr.R --expressions "$MET_EXPRS" --conditions "$MET_CONDITIONS" --state-1 "$STATE1" --state-2 "$STATE2" -o "$MET_PVALS"
#echo "Gene differential expression"
##./R/LogDiffExpr.R --expressions "$GENE_EXPRS" --conditions "$GENE_CONDITIONS" --deseq --state-1 "$STATE1" --state-2 "$STATE2" -o "$GENE_PVALS"
#./R/LogDiffExpr.R --expressions "$GENE_EXPRS" --conditions "$GENE_CONDITIONS" --log2 --state-1 "$STATE1" --state-2 "$STATE2" -o "$GENE_PVALS"

echo "Converting gene p-values to reaction p-values"
./R/gene2rxn_pval.R \
    --rxn2genes "$NETWORK_BASE.rxn2genes.tsv" \
    --gene-pvals "$GENE_PVALS" \
    -o "$RXN_PVALS"

tail -n +2 "$MET_PVALS" | \
cat \
    "$RXN_PVALS" \
    - \
    > "$WORKDIR/combined.pval"


#echo "Extracting subnet specific to p-values"
#./R/subnet.R \
#    -i "$WORKDIR/combined.pval" \
#    -n "$NETWORK_BASE.squared.sif" \
#    -o "$WORKDIR/combined.pval.specific.tab"

#echo "Analysing whole network"
#./network_analyis.R \
#    --fdr "$FDR" \
#    -i "$WORKDIR/combined.noncomp.pval" \
#    -n "$NETWORK_BASE.squared.noncomp.tab" \
#    -o "$WORKDIR/combined.noncomp.pval.module"

echo "Analysing noEx network"
./R/network_analyis.R \
    --fdr "$FDR" \
    -i "$WORKDIR/combined.pval" \
    -n "$NETWORK_BASE.nogene.masked.squared.nocomp.noex.hmdb.esc" \
    -o "$WORKDIR/combined.noncomp.pval.noEx.module" \
    --heinz ./heinz.py


#echo "Analysing metabolite network"
#./network_analyis.R \
#    --fdr "$FDR" \
#    -i "$MET_NONCOMP_PVALS" \
#    -n "$NETWORK_BASE.squared.noncomp.noEx.tab" \
#    -o "$WORKDIR/mets.noncomp.pval.noEx.module"
#
#echo "Analysing gene network"
#./network_analyis.R \
#    --fdr "$FDR" \
#    -i "$GENE_PVALS" \
#    -n $NETWORK_BASE.squared.noncomp.noEx.tab \
#    -o "$WORKDIR/genes.pval.noEx.module"
#
