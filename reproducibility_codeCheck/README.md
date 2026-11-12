# Container code checks

This directory contains small code execution checks for the MBv R and Python container environments. 
They are separate from the image build and smoke tests: these scripts exercise project code 
with public or project data and write analysis outputs. These are functional checks on reduced 
inputs. They do not validate full-data results or replace the pipeline's full analysis. 

These checks occurred prior to Nextflow implementation. These checks have not yet been designed for
portable image builds in a clean Linux setup.

<br>

## Python AUCell checks

The Python container check runs AUCell on a small subset of loom input. The full loom has 
21,265 cells. A full-data AUCell run was killed by Docker Desktop's 8 GiB memory limit 
(`OOMKilled=true`, exit status 137), even after reducing the worker count. Therefore, the 
Python scripts take the first 10% of cells (2,126 cells) and run with one worker.

### Python inputs and locations

The AUCell script requires access to a loom expression file and the `module.dat` file. 
The source files are found here: 
`/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_DEG_GRN/`
The Docker check requires these files to be staged locally while the Apptainer check 
from JHPCE accesses the originals on `/dcs04/`. These small differences account for the 
separate `.py` scripts for each container build.

Prior to the running the Docker check, SFTP transfer these files into the `reproducibility_codeCheck` 
under a new `data` dir, mirroring the original directory structure containing an `expr_loom` subdir: 

- `reproducibility_codeCheck/data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat`
- `reproducibility_codeCheck/data/expr_loom/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom`

The files are not committed; `.gitignore` excludes `reproducibility_codeCheck/data/`.

### Run the Docker Python check

Requirements: Docker Desktop (or Docker Engine), the built `mbv-python:test` image for
`linux/amd64`, the two local inputs above, and a project-root working directory.

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

The run writes a log to `test_results/docker/aucell.log` and a successful run creates an AUCell
loom file in the same location. 

### Run the JHPCE Apptainer Python check

Requirements: JHPCE access, `/dcs04/` access, the Python JHPCE SIF **at the path set in 
`apptainer_python.sh`**, and Apptainer (loaded through the pre-built module). 

Submit from the project root so Slurm resolves its relative output path correctly:

```bash
mkdir -p reproducibility_codeCheck/data/expr_loom test_results/apptainer
sbatch reproducibility_codeCheck/apptainer_python.sh
```

The run writes a Slurm log to `test_results/apptainer/mbv-aucell-python_<job-id>.log` and 
a successful run creates an AUCell loom file in the same location. 

<br><br>

## R GEOquery & SPE build checks

The R container check downloads public GEO data from the MBv accession `GSE333974` into a temporary directory and builds a small `SpatialExperiment` (SPE) from three samples. The same `.R` script is used for both Docker and Apptainer checks.

`Rscript --vanilla` skips the mounted checkout's `.Rprofile`, which would activate its local renv library instead of using the package library installed in the image.

### Run the Docker R check

Requirements: Docker Desktop (or Docker Engine), the built `mbv-r:test` image, network access to NCBI GEO, and enough Docker Desktop storage for the archive and its extracted contents.

Run from the MBv project root in zsh:

```zsh
mkdir -p test_results/docker
set -o pipefail
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" \
  -w /work \
  mbv-r:test \
  Rscript --vanilla reproducibility_codeCheck/spe.R test_results/docker/ \
  2>&1 | tee test_results/docker/spe.log
print "exit status: $?"
```
The project bind preserves the SPE output and log on the host; `--rm` discards the container.
The run writes a log to `test_results/docker/spe.log` and a successful run creates an HDF5 SPE
representation that consists of two, linked files in the same location. 

### Run the Apptainer R check

Requirements: JHPCE access, the R JHPCE SIF **at the path set in `apptainer_r.sh`**, and 
Apptainer (loaded through the pre-built module). 

Submit from the project root, with the Slurm log directory present before submission:

```bash
mkdir -p test_results/apptainer
sbatch reproducibility_codeCheck/apptainer_r.sh
```

The run writes a Slurm log to `test_results/apptainer/mbv-spe-r_<job-id>.log` and a successful 
run creates an HDF5 SPE representation that consists of two, linked files in the same location.

<br><br>

## Outputs at a glance

| Output | Docker Python | JHPCE Apptainer Python | Docker R | JHPCE Apptainer R |
| --- | --- | --- | --- | --- |
| Run log | `test_results/docker/aucell.log` | `test_results/apptainer/mbv-aucell-python_<job-id>.log` | `test_results/docker/spe.log` | `test_results/apptainer/mbv-spe-r_<job-id>.log` |
| Temporary/subset input | OS temporary directory (`mbv-aucell-*/*_first10pct.loom`); removed when Python exits | OS temporary directory (`mbv-aucell-*/*_first10pct.loom`); removed when Python exits | GEO download/extraction under R `tempdir()`; removed at script end | GEO download/extraction under R `tempdir()`; removed at script end |
| Analysis output | `test_results/docker/*_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` | `test_results/apptainer/*_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` | `test_results/docker/spe_n3_se.rds` and `test_results/docker/spe_n3_assays.h5` | `test_results/apptainer/spe_n3_se.rds` and `test_results/apptainer/spe_n3_assays.h5` |
