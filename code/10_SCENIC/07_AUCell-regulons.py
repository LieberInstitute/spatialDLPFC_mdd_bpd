import os, sys, pickle

from shutil import copyfile
from pyscenic.aucell import aucell
from ctxcore.genesig import openfile

from pyscenic.cli.utils import (
    ATTRIBUTE_NAME_CELL_IDENTIFIER,
    ATTRIBUTE_NAME_GENE,
    append_auc_mtx,
    is_valid_suffix,
    load_adjacencies,
    load_exp_matrix,
    load_modules,
    load_signatures,
    save_enriched_motifs,
    save_matrix,
    suffixes_to_separator,
)

ATTRIBUTE_NAME_WORKERS = 10
ATTRIBUTE_NAME_SEED = 1234

# passed args
cluster = sys.argv[1]

SUBSET="seurat-label-"+cluster
SUBSET_ID="spe-n109_"+SUBSET+"_13162-genes_no-lowUMI"

pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
RESULTS_DIR = pdir+"processed-data/10_SCENIC/"
LOOM_DIR = RESULTS_DIR+"expr_loom/"
REG_DIR = RESULTS_DIR+"AUCell_regulons/"

REGULONS_DAT_FNAME = os.path.join(RESULTS_DIR, "spe-n109_13162-no-lowUMI_logcounts.regulons-top20.dat")
LOOM_FNAME = os.path.join(LOOM_DIR, '{}_logcounts.loom'.format(SUBSET_ID))
OUT_FNAME = os.path.join(REG_DIR, '{}_logcounts_regulons-top20_AUCell-fixed-thold-05.loom'.format(SUBSET_ID))


ex_mtx = load_exp_matrix(
            LOOM_FNAME,
            False, #transpose arg
            False,  # sparse loading is disabled here for now
            ATTRIBUTE_NAME_CELL_IDENTIFIER,
            ATTRIBUTE_NAME_GENE,
        )

## drop any genes where all values for subset of cells ==0
#countRows = (ex_mtx!=0).sum()
#min_rows = round(.01*len(ex_mtx.index.values))-1
#drop_mtx = ex_mtx.loc[:, countRows>min_rows]
#print("Dropped genes present in <1% of spots - new shape:", drop_mtx.shape)

##set auc threshold so that max rank is equal(ish) to 10% detected genes cutoff
#auc_thold = round(nGenes/len(drop_mtx.columns), 3)
#print("Max rank of", nGenes, "is approx equal to an AUC threshold of:", auc_thold)

#set auc threshold so that max rank is 450 which for 13844 genes is 0.033
#auc_thold =.033
auc_thold = .05
print("AUC threshold set to 0.05 for all clusters")

with open(REGULONS_DAT_FNAME, 'rb') as file:
    signatures = pickle.load(file)

auc_mtx = aucell(
	ex_mtx,
        signatures,
        auc_threshold=auc_thold,
        noweights=False, # must set to False so that adjacency weights (aka importance) are used in enrichment
        seed=ATTRIBUTE_NAME_SEED,
        num_workers=ATTRIBUTE_NAME_WORKERS,
    )

print(auc_mtx.shape)
print(auc_mtx.head())

copyfile(LOOM_FNAME, OUT_FNAME)
append_auc_mtx(OUT_FNAME,
                ex_mtx,
                auc_mtx,
                signatures,
                ATTRIBUTE_NAME_SEED,
                ATTRIBUTE_NAME_WORKERS,
            )
