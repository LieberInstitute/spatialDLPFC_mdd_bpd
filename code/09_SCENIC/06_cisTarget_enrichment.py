import os, sys, pickle

import logging
from pyscenic.log import create_logging_handler
from typing import Sequence, Type
import pandas as pd

from ctxcore.rnkdb import RankingDatabase, opendb
from ctxcore.genesig import openfile
from pyscenic.prune import _prepare_client, find_features, prune2df
from pyscenic.transform import df2regulons

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

LOGGER = logging.getLogger(__name__)

pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
RESULTS_DIR = pdir+"processed-data/09_SCENIC/"
#LOOM_DIR = RESULTS_DIR+"expr_loom/"
#MOD_DIR = RESULTS_DIR+"AUCell_modules/"


ANNOT_PATH=[pdir+"raw-data/SCENIC_aux/genome_annotation/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
	pdir+"raw-data/SCENIC_aux/genome_annotation/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
MOTIF_PATH=pdir+"raw-data/SCENIC_aux/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

DATASET_ID='spe-n119_13162-no-lowUMI_logcounts'


MODULES_DAT_FNAME = os.path.join(RESULTS_DIR, '{}.modules-top20.dat'.format(DATASET_ID))
REGULONS_FNAME = os.path.join(RESULTS_DIR, '{}_regulons-weighted-top20.csv'.format(DATASET_ID))
REGULONS_DAT_FNAME = os.path.join(RESULTS_DIR, '{}.regulons-weighted-top20.dat'.format(DATASET_ID))

with open(MODULES_DAT_FNAME, 'rb') as file:
#with open('processed-data/09_SCENIC/spe-n119_13162-no-lowUMI_logcounts.modules.dat', 'rb') as file:
    signatures = pickle.load(file)


LOGGER.info("Loading databases.")
#def _load_dbs(fnames: Sequence[str]) -> Sequence[Type[RankingDatabase]]:
#    def get_name(fname):
#        return os.path.splitext(os.path.basename(fname))[0]
#
#    return [opendb(fname=fname.name, name=get_name(fname.name)) for fname in fnames]

dbs = [opendb(fname=ANNOT_PATH[0], name="promoter"),
	opendb(fname=ANNOT_PATH[1], name="enhancer")]


calc_func = prune2df

df_motifs = calc_func(
  dbs,
  signatures,
  MOTIF_PATH,
  rank_threshold= 5000,
  auc_threshold= .05,
  nes_threshold= 3.0,
  client_or_address= "custom_multiprocessing",
  module_chunksize= 100,
  num_workers= 10,
  motif_similarity_fdr= 1e-5,
  orthologuous_identity_threshold= 0.8,
  weighted_recovery=True,
)
print("Unfiltered regulons:", df_motifs.shape)

# filter out regulons that aren't directly annotated or have qvalue==0
df_motifs.columns = df_motifs.columns.droplevel(0)
keep1 = df_motifs['Annotation'].str.contains('directly annotated')
keep2 = df_motifs['MotifSimilarityQvalue']==0
df_motifs = df_motifs[keep1 | keep2]

print("Regulons (filtered for direct annotation or q==0):", df_motifs.shape)

LOGGER.info("Writing results to csv file.")
#df_motifs.to_csv(RESULTS_DIR+DATASET_ID+"_regulons-weighted.csv")
df_motifs.to_csv(REGULONS_FNAME)

LOGGER.info("Pickling regulons to .dat file.")
regulons = df2regulons(df_motifs)
with open(REGULONS_DAT_FNAME, 'wb') as f:
#with open("/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_SCENIC/spe-n119_13844-no-lowUMI_logcounts.regulons-weighted.dat", 'wb') as f:
    pickle.dump(regulons, f)
