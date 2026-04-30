# Python venv setup for tensorQTL (`uv`)

This documents how to recreate the local `.venv` for this project from the
pinned package set in `requirements.txt`.

All commands assume the current directory is `11_eQTL_coloc`.

## 1) Prerequisites

- `uv` installed (`uv 0.9.28` used for this environment)
- network access to:
  - PyPI
  - PyTorch CUDA index: `https://download.pytorch.org/whl/cu126`

## 2) Refresh pinned requirements from the current `.venv`

If the current `.venv` has been intentionally updated, refresh the tracked
requirements file before recreating the environment elsewhere:

```bash
uv pip freeze --python .venv/bin/python | sort > requirements.txt
```

The current file contains 98 pinned packages.

## 3) Create a local uv venv

```bash
uv venv --python 3.12.8 .venv
```

## 4) Install exact pinned requirements into `.venv`

`torch==2.9.1+cu126` and `torchvision==0.24.1+cu126` require the PyTorch CUDA
wheel index. `uv` must be allowed to resolve across both indexes:

```bash
uv pip sync --python .venv/bin/python requirements.txt \
  --extra-index-url https://download.pytorch.org/whl/cu126 \
  --index-strategy unsafe-best-match
```

## 5) Apply tensorQTL patch

Patch file in this project:

- `tensorqtl_core_invex.patch`

This patch fixes a minor tensorQTL failure mode in `tensorqtl/core.py` during
interaction eQTL testing. The original code uses batched `.inverse()` on the
per-variant regression matrix. If one variant produces a singular design
matrix, that inversion can raise an error and stop the entire run.

The patch replaces that inversion with
`torch.linalg.inv_ex(..., check_errors=False)`, checks the returned `info`
values, and marks singular per-variant batch elements with `NaN` instead of
crashing the run.

Apply the patch inside the local site-packages root so the patch path
`tensorqtl/core.py` matches:

```bash
cd .venv/lib/python3.12/site-packages
patch -p0 < ../../../../tensorqtl_core_invex.patch
cd -
```

## 6) Verify environment and patch

```bash
.venv/bin/python - <<'PY'
from importlib.metadata import version
mods = ['tensorqtl', 'torch', 'pandas', 'rpy2', 'numpy', 'scipy', 'pyarrow']
for m in mods:
    print(m, version(m))
PY
```

Expected key versions:

- `tensorqtl==1.0.10`
- `torch==2.9.1+cu126`
- `pandas==2.3.3`
- `rpy2==3.6.4`

Patch verification:

```bash
sed -n '212,252p' .venv/lib/python3.12/site-packages/tensorqtl/core.py
```

You should see code using `XtX`, `torch.linalg.inv_ex(...)`, and singular-batch
handling via `info > 0`.
