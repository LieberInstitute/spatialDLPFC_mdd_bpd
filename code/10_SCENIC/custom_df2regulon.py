import pandas as pd

from ctxcore.genesig import Regulon
from typing import Optional, Sequence, Type
from ctxcore.rnkdb import RankingDatabase
from pyscenic.transform import _regulon4group

from pyscenic.utils import (
    ACTIVATING_MODULE,
    COLUMN_NAME_ANNOTATION,
    COLUMN_NAME_MOTIF_ID,
    COLUMN_NAME_MOTIF_SIMILARITY_QVALUE,
    COLUMN_NAME_ORTHOLOGOUS_IDENTITY,
    COLUMN_NAME_TF,
    REPRESSING_MODULE,
)

COLUMN_NAME_NES = "NES"
COLUMN_NAME_AUC = "AUC"
COLUMN_NAME_CONTEXT = "Context"
COLUMN_NAME_TARGET_GENES = "TargetGenes"
COLUMN_NAME_RANK_AT_MAX = "RankAtMax"
COLUMN_NAME_TYPE = "Type"

def df2regulons_custom(df, save_columns=[]) -> Sequence[Regulon]:
    df = df.copy()

    if df.columns.nlevels == 2:
        df.columns = df.columns.droplevel(0)

    def get_type(row):
        ctx = row[COLUMN_NAME_CONTEXT]
        return 'repressing' if 'enhancer' in ctx else 'activating'

    df[COLUMN_NAME_TYPE] = df.apply(get_type, axis=1)

    not_none = lambda r: r is not None
    return list(
        filter(
            not_none,
            (
                _regulon4group(
                    tf_name, frozenset([interaction_type]), df_grp, save_columns
                )
                for (tf_name, interaction_type), df_grp in df.groupby(
                    by=[COLUMN_NAME_TF, COLUMN_NAME_TYPE]
                )
            ),
        )
    )

