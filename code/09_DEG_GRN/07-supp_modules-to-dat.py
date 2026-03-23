import os, pickle
from functools import partial
from itertools import chain
from typing import Sequence, Type
from urllib.parse import urljoin
import attr

import numpy as np
import pandas as pd
from ctxcore.genesig import GeneSignature, Regulon, openfile

from pyscenic.utils import modules4thr
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

pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
RESULTS_DIR = pdir+"processed-data/09_DEG_GRN/"

# input file
ADJ_FNAME = RESULTS_DIR+"spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv"

# output files
MODULES_FNAME = RESULTS_DIR+"spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.csv"
MODULES_DAT_FNAME = RESULTS_DIR+"spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat"

# load adjacencies
adjacencies = load_adjacencies(ADJ_FNAME)
print("Adj. matrix dims:", adjacencies.shape)

def select_modules_from_adjacencies(
    adjacencies: pd.DataFrame,
    thresholds=(0.75, 0.90),
    min_genes=20,
    absolute_thresholds=True,
    keep_only_activating=False
) -> Sequence[Regulon]:

    if not absolute_thresholds:

        def iter_modules(adjc, context):
            yield from chain(
                chain.from_iterable(
                    modules4thr(
                        adjc, thr, context, pattern="weight>{}%".format(frac * 100)
                    )
                    for thr, frac in zip(
                        list(adjacencies["importance"].quantile(thresholds)),
                        thresholds,
                    )
                ),
            )

    else:

        def iter_modules(adjc, context):
            yield from chain(
                chain.from_iterable(
                    modules4thr(adjc, thr, context) for thr in thresholds
                ),
            )


    activating_modules = adjacencies[adjacencies['regulation'] > 0.0]
    if keep_only_activating:
        modules_iter = iter_modules(
            activating_modules, frozenset(['activating'])
        )
    else:
        repressing_modules = adjacencies[adjacencies['regulation'] < 0.0]
        modules_iter = chain(
            iter_modules(activating_modules, frozenset(['activating'])),
            iter_modules(repressing_modules, frozenset(['repressing'])),
        )

    # Derive modules for these adjacencies.
    # + Add the transcription factor to the module.
    #   [We are unable to assess if a TF works in a direct self-regulating way, either inhibiting its own expression or
    #    activating it. Therefore the most unbiased way forward is to add the TF to both activating as well as
    #    repressing modules]
    # + Filter for minimum number of genes.
    #LOGGER.info("Creating modules.")

    def add_tf(module):
        return module.add(module.transcription_factor)

    return list(filter(lambda m: len(m) >= min_genes, map(add_tf, modules_iter)))

modules = select_modules_from_adjacencies(adjacencies, thresholds=(1,),
	min_genes=10, absolute_thresholds=True,)

# save as DF for easy browsing
mod_df = pd.DataFrame({"TF": [item.transcription_factor for item in modules],
	"dir": ["activating" if "activating" in item.context else "repressing" for item in modules],
	"set_size": [len(item) for item in modules],
	"set_str": ['/' .join(item.genes) for item in modules]
})

print("Module df dims:", mod_df.shape)
mod_df.to_csv(MODULES_FNAME)

with open(MODULES_DAT_FNAME, 'wb') as f:
    pickle.dump(modules, f)
