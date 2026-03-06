import os, pickle

from shutil import copyfile
from pyscenic.aucell import aucell

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

# auc command
SUBSET="seurat-label-Micro.Vasc"
SUBSET_ID="spe-n119_"+SUBSET+"_13844-genes_no-lowUMI"

pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
RESULTS_DIR = pdir+"processed-data/09_SCENIC/"
LOOM_DIR = RESULTS_DIR+"expr_loom/"
MOD_DIR = RESULTS_DIR+"AUCell_modules/"

#MODULES_DAT_FNAME = os.path.join(RESULTS_DIR, "spe-n119_13844-no-lowUMI_logcounts.modules.dat")
LOOM_FNAME = os.path.join(LOOM_DIR, '{}_logcounts.loom'.format(SUBSET_ID))
OUT_FNAME = os.path.join(MOD_DIR, '{}_logcounts_modules-AUCell.loom'.format(SUBSET_ID))


ex_mtx = load_exp_matrix(
            LOOM_FNAME,
            False, #transpose arg
            False,  # sparse loading is disabled here for now
            ATTRIBUTE_NAME_CELL_IDENTIFIER,
            ATTRIBUTE_NAME_GENE,
        )

with open('processed-data/09_SCENIC/spe-n119_13844-no-lowUMI_logcounts.modules.dat', 'rb') as file:
    signatures = pickle.load(file)

auc_mtx = aucell(
        ex_mtx,
        signatures,
        auc_threshold=.05,
        noweights=False,
        seed=1234,
        num_workers=ATTRIBUTE_NAME_WORKERS,
    )

print(auc_mtx.head())

copyfile(LOOM_FNAME, OUT_FNAME)
append_auc_mtx(OUT_FNAME,
                ex_mtx,
                auc_mtx,
                signatures,
                1234,
                ATTRIBUTE_NAME_WORKERS,
            )
