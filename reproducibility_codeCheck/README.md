# Python AUCell code check

This directory contains reduced-input checks that run the project's AUCell analysis code in the MBv Python container. These are **code execution checks**, separate from the container build and smoke tests: they load the project inputs, calculate AUCell scores, and write an output loom file.

The same analysis is run in two environments:

- **Docker:** `docker-aucell.py` reads inputs staged under `reproducibility_codeCheck/data/` in the project checkout.
- **JHPCE Apptainer:** `apptainer-aucell.py` reads the source loom and module data from the JHPCE project data directory, mounted read-only at `/source-data`.

Both tests write the reduced loom to the project checkout and keep their AUCell result and run log in separate Docker and Apptainer result folders.

## Why the test uses the first 10% of cells

The full loom has 21,265 cells. A full-data AUCell run was killed by Docker Desktop's 8 GiB memory limit (`OOMKilled=true`, exit status 137), even after reducing the worker count. The test scripts therefore take the first 10% of cells (2,126 cells) and run with one worker to make a lightweight functional check.

This subset checks that the container can read the loom and module data and execute AUCell. It is not a full-data pipeline run, a measurement of full-data memory requirements, or a validation of scientific results on all cells. The scripts create a separate `*_first10pct.loom` file; they do not alter the source loom. Loom compression means the subset file size need not be exactly 10% of the original.

The check uses the `Micro.Vasc` cluster's module signatures, AUCell threshold `0.05`, seed `1234`, and weighted scoring (`noweights=False`).

## Inputs and their locations

### Docker

Run from the MBv project root. The following inputs must be available locally:

- `reproducibility_codeCheck/data/expr_loom/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom`
- `reproducibility_codeCheck/data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat`

These input data files are not committed: `.gitignore` excludes `data/`. Make sure the files are available at the paths above before running the Docker check.

### JHPCE Apptainer

The batch script `apptainer.sh` expects the source files under:

`/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_DEG_GRN/`

It binds only that directory, read-only, to `/source-data` in the container. The expected files are:

- `/source-data/expr_loom/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom`
- `/source-data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat`

The project checkout is bound to `/work`, and the Apptainer script writes the first-10%-cells loom under `/work/reproducibility_codeCheck/data/expr_loom/`.

## Docker check

Requirements: Docker Desktop (or Docker Engine), the built `mbv-python:test` image for `linux/amd64`, the two local inputs above, and `test_results/docker/`.

From the project root, run in zsh:

```zsh
# Add Docker's CLI to PATH if this shell cannot find `docker`.
export PATH="/Applications/Docker.app/Contents/Resources/bin:$PATH"
docker info

# Capture output and preserve the Docker command's exit status through tee.
mkdir -p test_results/docker
set -o pipefail
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" \
  -w /work \
  mbv-python:test \
  python3 -u reproducibility_codeCheck/docker-aucell.py \
  2>&1 | tee test_results/docker/aucell.log
print "exit status: $?"
```

The bind mount maps the checkout to `/work`; generated files therefore persist in the checkout. `tee` replaces `test_results/docker/aucell.log` on each run.

## JHPCE Apptainer check

Requirements: JHPCE access, the built Python JHPCE SIF at the path set in `apptainer.sh`, the two source files under `processed-data/09_DEG_GRN/`, and a compute allocation with Apptainer available. Submit from the MBv project root so the relative Slurm output path resolves correctly.

```bash
mkdir -p reproducibility_codeCheck/data/expr_loom test_results/apptainer
sbatch reproducibility_codeCheck/apptainer.sh
```

The script loads `apptainer/1.5.3`, binds the project checkout to `/work`, and binds only `processed-data/09_DEG_GRN` to `/source-data` as read-only. Its `#SBATCH --output` file captures both standard output and errors, including the final exit status; it does not use a second `tee` log.

## Outputs and interpreting results

| Output | Docker | JHPCE Apptainer |
| --- | --- | --- |
| Run log | `test_results/docker/aucell.log` | `test_results/apptainer/mbv-aucell-python_<job-id>.log` |
| First-10%-cells loom | `reproducibility_codeCheck/data/expr_loom/*_first10pct.loom` | `reproducibility_codeCheck/data/expr_loom/*_first10pct.loom` |
| AUCell result loom | `test_results/docker/*_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` | `test_results/apptainer/*_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` |

Exit status `0` means the script completed successfully. Inspect the corresponding log for the subset cell count, AUCell output shape and example scores, and final exit status. A nonzero status indicates failure; check the traceback or error in the log. For Docker, status `137` together with `OOMKilled=true` indicates the container exceeded its memory limit.

The older `reproducibility_codeCheck/aucell-docker.log` is a saved earlier Docker capture; the command above writes the current run to `test_results/docker/aucell.log`.

These checks cover the Python image and a reduced AUCell run only. They do not test the R image, portable Apptainer build, full 21,265-cell AUCell run, or full pipeline outputs.
