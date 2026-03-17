import os, pickle

import pandas as pd
import ast

from ctxcore.genesig import openfile
from pyscenic.transform import df2regulons

pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
RESULTS_DIR = pdir+"processed-data/10_SCENIC/"

DATASET_ID='spe-n119_13162-no-lowUMI_logcounts'
REGULONS_FNAME = os.path.join(RESULTS_DIR, '{}_regulons-weighted-top20_refined-filtered.csv'.format(DATASET_ID))
REGULONS_DAT_FNAME = os.path.join(RESULTS_DIR, '{}.regulons-weighted-top20-refined-filtered.dat'.format(DATASET_ID))

df_motifs = pd.read_csv(REGULONS_FNAME)
print("\nRefined regulons before merge:", df_motifs.shape)

def str_to_list(row):
    tgs = row['TargetGenes']
    return ast.literal_eval(tgs)

df_motifs['TargetGenes'] = df_motifs.apply(str_to_list, axis=1)

regulons = df2regulons(df_motifs)
print("\nRefined regulons afer merge:", len(regulons))

with open(REGULONS_DAT_FNAME, 'wb') as f:
    pickle.dump(regulons, f)
