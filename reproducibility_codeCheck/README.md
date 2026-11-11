# Container code checks

This directory contains small code execution checks for the MBv R and Python container environments. They are separate from the image build and smoke tests: these scripts exercise project code with public or project data and write analysis outputs.

The checks currently documented here are:

- **Docker Python:** runs AUCell on the first 10% of cells in a loom input.
- **JHPCE Apptainer Python:** runs the same AUCell check using source data from JHPCE.
- **Docker R:** downloads public GEO data and builds a small `SpatialExperiment` (SPE) from three samples.

These are functional checks on reduced inputs. They do not validate full-data results or replace the pipeline's full analysis.

## Python AUCell checks

### Why the Python test uses 10% of cells

The full loom has 21,265 cells. A full-data AUCell run was killed by Docker Desktop's 8 GiB memory limit (`OOMKilled=true`, exit status 137), even after reducing the worker count. The Python scripts therefore take the first 10% of cells (2,126 cells) and run with one worker.

This subset checks that the container can read the loom and module data and execute AUCell. It is not a full-data run, a measurement of full-data memory requirements, or a validation of scientific results on all cells. The scripts create a separate `*_first10pct.loom` file and leave the source loom unchanged. Loom compression means the subset file size need not be exactly 10% of the original.

The AUCell check uses the `Micro.Vasc` cluster's module signatures, threshold `0.05`, seed `1234`, and weighted scoring (`noweights=False`).

### Python inputs and locations

#### Docker Python

Run from the MBv project root. These input data files must be staged locally:

- `reproducibility_codeCheck/data/expr_loom/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom`
- `reproducibility_codeCheck/data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat`

The files are not committed; `.gitignore` excludes `reproducibility_codeCheck/data/`.

#### JHPCE Apptainer Python

The batch script `apptainer.sh` reads these source files from:

`/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_DEG_GRN/`

It binds only this directory, read-only, to `/source-data` in the container. The expected files are:

- `/source-data/expr_loom/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom`
- `/source-data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat`

The checkout is bound to `/work`; the first-10%-cells loom is written under `/work/reproducibility_codeCheck/data/expr_loom/`.

### Run the Docker Python check

Requirements: Docker Desktop (or Docker Engine), the built `mbv-python:test` image for `linux/amd64`, the two local inputs above, and a project-root working directory.

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

The checkout bind makes inputs and outputs available on the host. `tee` replaces `test_results/docker/aucell.log` on each run.

### Run the JHPCE Apptainer Python check

Requirements: JHPCE access, the Python JHPCE SIF at the path set in `apptainer.sh`, the two source files under `processed-data/09_DEG_GRN/`, and a compute allocation with Apptainer available. Submit from the project root so Slurm resolves its relative output path correctly:

```bash
mkdir -p reproducibility_codeCheck/data/expr_loom test_results/apptainer
sbatch reproducibility_codeCheck/apptainer.sh
```

The script loads `apptainer/1.5.3`, binds the checkout to `/work`, and binds only `processed-data/09_DEG_GRN` to `/source-data` as read-only. The `#SBATCH --output` file captures standard output and errors, including the exit status; it does not create a second `tee` log.

## Docker R GEOquery and SPE check

This code check exercises the R image's GEOquery, SpatialExperiment, and supporting packages while retrieving public data and creating an SPE from the first three samples. It runs `reproducibility_codeCheck/docker-spe.R`.

### Data source and temporary files

The script retrieves dataset **GSE333974** from NCBI GEO using `GEOquery::getGEO()` and `getGEOSuppFiles()`. No local copy of the GEO input files is required, but the container needs network access to NCBI GEO.

The script downloads the full `GSE333974_RAW.tar` archive (about 4.2 GB in the recorded run), extracts its contents into an R `tempdir()`, and uses the first three samples to construct the SPE. The temporary directory is inside the container and is removed by the script at the end; it does not need its own bind mount. The full archive and extracted contents require temporary disk space during the run, even though only three samples are used for the SPE.

The persistent outputs are written under `test_results/docker/`, which is part of the bound project checkout.

### Run the Docker R check

Requirements: Docker Desktop (or Docker Engine), the rebuilt `mbv-r:test` image containing the updated R lockfile, network access to NCBI GEO, and enough Docker Desktop storage for the archive and its extracted contents.

Run from the MBv project root in zsh:

```zsh
mkdir -p test_results/docker
set -o pipefail
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" \
  -w /work \
  mbv-r:test \
  Rscript --vanilla reproducibility_codeCheck/docker-spe.R \
  2>&1 | tee test_results/docker/spe.log
print "exit status: $?"
```

`Rscript --vanilla` skips the mounted checkout's `.Rprofile`, which would activate its local renv library instead of using the package library installed in the image. The project bind preserves the SPE output and log on the host; `--rm` discards the container and its temporary download files after exit. `tee` replaces `test_results/docker/spe.log` on each run.

### R outputs and pass/fail

- `test_results/docker/spe.log` records GEOquery downloads, SPE construction, timing, session information, and completion marker `*** TEST END ***`.
- `test_results/docker/spe_n3_se.rds` and `test_results/docker/spe_n3_assays.h5` are the HDF5-backed SPE outputs produced with prefix `spe_n3_`.
- Exit status `0` indicates the script completed. A nonzero status indicates failure; inspect the log for the error. The log includes separate timing for `read10xVisium()` and saving the HDF5-backed object, but it does not record total script runtime.

## Outputs at a glance

| Output | Docker Python | JHPCE Apptainer Python | Docker R |
| --- | --- | --- | --- |
| Run log | `test_results/docker/aucell.log` | `test_results/apptainer/mbv-aucell-python_<job-id>.log` | `test_results/docker/spe.log` |
| Temporary/subset input | `reproducibility_codeCheck/data/expr_loom/*_first10pct.loom` | `reproducibility_codeCheck/data/expr_loom/*_first10pct.loom` | GEO download/extraction under R `tempdir()`; removed at script end |
| Analysis output | `test_results/docker/*_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` | `test_results/apptainer/*_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` | `test_results/docker/spe_n3_se.rds` and `test_results/docker/spe_n3_assays.h5` |

The older `reproducibility_codeCheck/aucell-docker.log` is a saved Docker Python capture; the current Python command writes to `test_results/docker/aucell.log`.

These checks cover the Python image with a reduced AUCell input and the Docker R image with public GEO data and three samples. They do not test the R Apptainer code path, portable Apptainer builds on a clean Linux host, a full 21,265-cell AUCell run, or full-scale SPE construction from every GEO sample.
