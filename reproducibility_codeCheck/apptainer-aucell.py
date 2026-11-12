# first subset loom to 10% to try and fix docker memory issue
import atexit
import tempfile
from pathlib import Path
import loompy

#pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
#RESULTS_DIR = "processed-data/09_DEG_GRN/"
# pdir + RESULTS_DIR is now saved a /source-data/ within the container
source_loom = Path(
    "/source-data/expr_loom/"
    "spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom"
)
subset_tmp = tempfile.TemporaryDirectory(prefix="mbv-aucell-")
atexit.register(subset_tmp.cleanup)
subset_dir = Path(subset_tmp.name)
subset_loom = subset_dir / (
	f"{source_loom.stem}_first10pct{source_loom.suffix}"
)

with loompy.connect(str(source_loom), mode="r") as ds:
    # Loom convention: genes are rows and cells are columns.
    n_cells = ds.shape[1]
    n_keep = max(1, n_cells // 10)
    cell_indices = slice(0, n_keep)

    loompy.create(
        str(subset_loom),
        ds[:, cell_indices],
        ds.ra,
        {key: values[cell_indices] for key, values in ds.ca.items()},
    )

print(f"Created temporary subset {subset_loom} with the first {n_keep} of {n_cells} cells.")



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

ATTRIBUTE_NAME_WORKERS = 1
ATTRIBUTE_NAME_SEED = 1234

# passed args
cluster = "Micro.Vasc"

SUBSET="seurat-label-"+cluster
SUBSET_ID="spe-n119_"+SUBSET+"_13162-genes_no-lowUMI_logcounts"

#pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
#RESULTS_DIR = "processed-data/09_DEG_GRN"
MOD_DIR = "test_results/apptainer/"

MODULES_DAT_FNAME = Path("/source-data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat")
LOOM_FNAME = str(subset_loom)
OUT_FNAME = os.path.join(MOD_DIR, '{}_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom'.format(SUBSET_ID))


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

with open(MODULES_DAT_FNAME, 'rb') as file:
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
