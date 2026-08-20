#!/usr/bin/env python3

"""Write the exact tensorQTL covariate QR bases for targeted LD."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--datasets", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for dataset_id in args.datasets.split(","):
        cov_path = input_dir / f"{dataset_id}.gene.covars.txt"
        cov = pd.read_csv(cov_path, sep="\t", index_col=0).T

        ## reproduce tensorqtl.core.Residualizer 1.0.10 exactly
        cov_t = torch.tensor(cov.values, dtype=torch.float32)
        centered = cov_t - cov_t.mean(0)
        q_t, r_t = torch.linalg.qr(centered)

        q_path = output_dir / f"{dataset_id}.tensorqtl_Q.tsv.gz"
        q = pd.DataFrame(q_t.numpy(), index=cov.index)
        q.to_csv(q_path, sep="\t", compression="gzip", index_label="sample_id")

        diag = torch.diag(r_t).detach().cpu().numpy()
        rows.append(
            {
                "dataset_id": dataset_id,
                "n_samples": cov.shape[0],
                "n_covariate_columns": cov.shape[1],
                "centered_matrix_rank": int(torch.linalg.matrix_rank(centered)),
                "q_columns": q_t.shape[1],
                "zero_centered_columns": int(
                    np.sum(np.linalg.norm(centered.numpy(), axis=0) == 0)
                ),
                "min_abs_r_diagonal": float(np.min(np.abs(diag))),
                "torch_version": torch.__version__,
                "q_file": str(q_path.resolve()),
            }
        )

    pd.DataFrame(rows).to_csv(
        output_dir / "tensorqtl_q_provenance.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
