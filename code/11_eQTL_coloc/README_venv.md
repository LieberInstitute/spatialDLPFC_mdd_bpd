# Python venv setup for tensorQTL (`uv`)

This documents the exact steps used to recreate a local `.venv` in this project using the same pinned package set as:

- `~/work/R/eqtl_MDD_OUD/.venv`

## 1) Prerequisites

- `uv` installed (`uv 0.9.28` used here)
- network access to:
  - PyPI
  - PyTorch CUDA index: `https://download.pytorch.org/whl/cu126`

## 2) Export exact package requirements from reference env

From `code/11_eQTL_coloc`:

```bash
uv pip freeze --python /home/gpertea/work/R/eqtl_MDD_OUD/.venv/bin/python \
  | sort > requirements_uv_eqtl_MDD_OUD.freeze.txt
```

This produced a pinned 98-package freeze file:

- `requirements_uv_eqtl_MDD_OUD.freeze.txt`

## 3) Create local uv venv with matching Python version

```bash
uv venv --python 3.12.8 .venv
```

## 4) Install exact pinned requirements into local `.venv`

`torch==2.9.1+cu126` and `torchvision==0.24.1+cu126` require the PyTorch CUDA wheel index, and uv must be allowed to resolve across both indexes:

```bash
uv pip sync --python .venv/bin/python requirements_uv_eqtl_MDD_OUD.freeze.txt \
  --extra-index-url https://download.pytorch.org/whl/cu126 \
  --index-strategy unsafe-best-match
```

## 5) Apply tensorQTL patch in local `.venv`

Patch file in this project:

- `tensorqtl_core_invex.patch`

Apply it inside local site-packages root (so patch path `tensorqtl/core.py` matches):

```bash
cd .venv/lib/python3.12/site-packages
patch -p0 < /home/gpertea/work/R/spatialDLPFC_mdd_bpd/code/11_eQTL_coloc/tensorqtl_core_invex.patch
```

## 6) Verify environment and patch

```bash
.venv/bin/python - <<'PY'
from importlib.metadata import version
mods=['tensorqtl','torch','pandas','rpy2','numpy','scipy','pyarrow']
for m in mods:
    print(m, version(m))
PY
```

Expected key versions (matching reference env):

- `tensorqtl==1.0.10`
- `torch==2.9.1+cu126`
- `pandas==2.3.3`
- `rpy2==3.6.4`

Patch verification (look for `torch.linalg.inv_ex` in `tensorqtl/core.py`):

```bash
sed -n '212,252p' .venv/lib/python3.12/site-packages/tensorqtl/core.py
```

You should see code using `XtX`, `torch.linalg.inv_ex(...)`, and singular-batch handling via `info > 0`.
