# MBv pipeline containers

To support consistent environment setup for analysis of the MBv spatial 
transcriptomics (SRT) data, and future application of the pipeline to external 
datasets, we have created the container images. These images support the bulk
of the current pipeline in code directories `00` through `10`. The tensorQTL work
in `11_eQTL_coloc` is outside the current scope and may receive a separate 
image later.

The MBv pipeline uses R and Python in separate jobs, so the environments are
provided as two images to keep their dependency sets isolated:

- **R image:** R 4.4.3 and the package versions recorded in `renv.lock`.
- **Python image:** the locked `pyscenic_bioconda` environment, Python 3.10.4,
  pySCENIC, and Arboreto. The pySCENIC install preserves a Dask memory patch
  required for full functionality with a large dataset.

<br>

Jump to [build instructions](#mbv-container-build-instructions).

<br>

## Container overview

Containers package the software environments used by the MBv pipeline so that
R and Python package versions, system libraries, and supporting tools can be
installed and run consistently across machines. They also make it easier for
users to run the pipeline without manually installing the original packages
one-by-one. 

The MBv pipeline containers recreate the recorded R and Python package environments
from the working setup on JHPCE (Joint High Performance Computing Exchange, the 
JHU Biostatistics Department's high-performance computing facility).

<br>

### Docker and Apptainer

Both container tools run Linux software images, but they serve different settings.
**Docker** is the convenient option for workstations, including macOS and
Windows through Docker Desktop. Desktop runs Linux containers in a managed
Linux VM. A user builds the images locally with Docker Desktop, and those images
are not automatically published to a registry. The MBv project images target 
`linux/amd64` to match the saved Python lock; Apple Silicon therefore uses 
emulation. 

**Apptainer** is intended for Linux systems and shared high-performance
computing environments. Its SIF images work with the site's scheduler and
filesystem setup, including JHPCE's Slurm environment. Apptainer definitions
are kept separate from Dockerfiles because the build tools and runtime model
differ.

<br>

### Portable and JHPCE Apptainer definitions

In building the Apptainer images within JHPCE, several site-specific workarounds 
were implemented. To improve application of the Apptainer definitions to other
high-performance computing environments, we created a portable Linux definition
and a JHPCE-specific definition:

| Environment | Portable Linux | JHPCE |
| --- | --- | --- |
| R | `apptainer/r-portable.def` | `apptainer/r-jhpce.def` |
| Python | `apptainer/python-portable.def` | `apptainer/python-jhpce.def` |

The portable definitions (yet to be tested) use normal Linux build behavior.
The JHPCE definitions apply site-specific workarounds required by its 
rootless/fakeroot build environment: APT runs with its sandbox user set to root, 
archive extraction uses `TAR_OPTIONS=--no-same-owner`, and staged verifier and R
lockfile permissions are normalized for runtime checks. Both profiles call
shared build scripts, so package installation logic stays aligned while the
JHPCE adjustments remain explicit in the JHPCE profile.

<br><br>

# MBv container build instructions

This guide explains how to build and smoke-test the separate R and Python
images for the bulk MBv pipeline environment in code directories `00`–`10`.
It covers container setup only, not implementation or execution of the
analysis workflow.

## Tested environments and availability

The JHPCE Apptainer and MacOS Docker R and Python images were built and smoke-tested 
successfully on 2026-11-09 (see details [below](#measured-image-sizes-times-and-memory)). 
R builds consistently take over an hour, so plan accordingly.

Docker images built with the commands below are local to that Docker
installation; they have not been published to a container registry.

<br>

Use this table to choose the build instructions best-suited for your setup.

| Machine/setup | Software needed to build and smoke-test | Project clone | Compatible image and current test status |
| --- | --- | --- | --- |
| Mac with Apple Silicon | Docker Desktop running; Terminal/Bash | Yes; clone `revision` | Docker `linux/amd64`; built and smoke-tested with emulation |
| Mac with Intel | Docker Desktop for Mac; Terminal/Bash | Yes; clone `revision` | Docker `linux/amd64`; matches CPU architecture, but this setup has not been tested |
| Windows PC with x86_64 CPU | Docker Desktop using Linux containers; WSL 2 and Bash terminal recommended | Yes; clone inside WSL filesystem | Docker `linux/amd64`; matches CPU architecture, but this setup has not been tested |
| Linux x86_64 with Docker | Docker Engine, Git, Bash | Yes; clone `revision` | Docker `linux/amd64`; clean Linux Docker build has not been tested |
| Linux x86_64 with Apptainer | Apptainer, Git, Bash; configured fakeroot or root access | Yes; clone `revision` | Portable Apptainer definitions are prepared; clean Linux builds have not been tested |
| JHPCE | JHPCE access, Slurm, and Apptainer module `apptainer/1.5.3` | Yes; checkout on `revision` | Separate named JHPCE Apptainer definitions; both built and passed runtime checks |

The Python lock targets Linux x86_64, so all Docker builds use
`linux/amd64`. Docker Desktop on Apple Silicon emulates that architecture;
native `linux/arm64` images would need a separate Conda lock and validation.

<br>

## Build and test steps

Briefly, you need to clone the project and then use bash to build the container image(s)
and test the builds. These tests verify the environment, not the full pipeline.
The build logs and test results are stored in `reproducibility_outputs`
and described in [Reading build and test results](#reading-build-and-test-results).

The R smoke test compares every installed package with `renv.lock` and loads
representative packages such as sf, s2, Seurat, SingleCellExperiment, and
spatialLIBD. The Python smoke test checks Python, pySCENIC, and Arboreto
versions, imports representative packages, and runs `pip check`. 

<br>

### Clone the project

Clone the `revision` branch of the MBv repo. The lockfiles, image definitions, 
and build scripts are all in this checkout; use it to build both images.

```bash
git clone --branch revision https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd.git
cd spatialDLPFC_mdd_bpd
```

#### Modified clone instructions

To omit cloning large files for just this container test, use these modified instructions
for a blobless partial clone with cone-mode sparse checkout. This keeps the project files
needed for environment setup and tests while leaving the top-level `raw-data` and 
`processed-data` directories out of the working tree:

```bash
git clone --filter=blob:none --sparse --branch revision \
  https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd.git \
  spatialDLPFC_mdd_bpd

cd spatialDLPFC_mdd_bpd

git ls-tree -d --name-only HEAD |
  grep -Ev '^(raw-data|processed-data)$' |
  git sparse-checkout set --stdin
```

Check the selection with:

```bash
git sparse-checkout list
test ! -e raw-data && echo "raw-data omitted"
test ! -e processed-data && echo "processed-data omitted"
git status --short
```
This only limits files materialized in this local checkout; it does not delete files from 
GitHub or alter branch history. It’s suitable for building and smoke-testing the images, 
but the MBv analysis itself requires the data.

<br>

### Build Option 1: Docker

Start Docker Desktop (or Docker Engine on Linux), then run the profiles
separately:

```bash
bash container/scripts/run_docker_test.sh r
bash container/scripts/run_docker_test.sh python
```

The runner builds each image, checks its locked versions and representative
imports, and repeats the checks in a fresh container. Docker images are tagged
locally as `mbv-r:test` and `mbv-python:test`. The runner saves commands, logs,
build and test durations, Docker versions, image ID, and image size under
`reproducibility_outputs/docker/<UTC-run-id>/<profile>/`.

On macOS, Docker Desktop must be running and have enough VM memory and disk.
On Windows, use Linux containers and run the commands from WSL 2.

<br>

### Build Option 2: Apptainer images on Linux

**These portable definitions have not yet been tested on a clean Linux host.**

Use a Linux x86_64 host with Apptainer installed, then run R and Python as
separate builds:

```bash
bash container/scripts/run_portable_test.sh r
bash container/scripts/run_portable_test.sh python
```

The runner defaults to `apptainer build --fakeroot`, which requires fakeroot
support configured in Apptainer. If the host has root access but no fakeroot
support, run the commands as root with `APPTAINER_BUILD_MODE=root`. The
portable definitions use normal Linux APT behavior and do not include the
JHPCE-specific settings. Logs, build and test durations, checksums, and
metadata go under `reproducibility_outputs/apptainer/<UTC-run-id>/<profile>-portable/`.

<br>

### Build Option 3: Apptainer images on JHPCE

The JHPCE definitions are `container/apptainer/r-jhpce.def` and
`container/apptainer/python-jhpce.def`. Use the login node to request a
separate 10 GiB compute allocation for each image. The runner loads the
Apptainer module, builds the selected image, and performs its smoke test.

```bash
ssh -X USERID@jhpce01.jhsph.edu
cd /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd
git switch revision
git pull origin revision

# First, build and test R in its own compute job.
srun --mem 10G --pty bash
cd /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd
export MBV_JHPCE_MEM_REQUEST=10G
bash container/scripts/run_jhpce_apptainer_test.sh r
exit

# Request a new compute allocation before building Python.
srun --mem 10G --pty bash
cd /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd
export MBV_JHPCE_MEM_REQUEST=10G
bash container/scripts/run_jhpce_apptainer_test.sh python
exit
exit
```

The JHPCE build profile uses the site's required APT sandbox setting,
`TAR_OPTIONS=--no-same-owner` during archive extraction, and readable
permissions for staged verifier files. The R profile also sets permissions on
`/opt/mbv/project` and its copied `renv.lock` so the unprivileged runtime test
can read the lockfile. These are the JHPCE-specific adjustments; the portable
profile does not use them.

The runner writes the selected definition, commands, logs, image checksum,
size, and build/test timing under
`reproducibility_outputs/apptainer/<UTC-run-id>/<profile>-jhpce/`. After each
compute job exits, capture peak memory from the login node and save the Slurm
record in that profile directory:

```bash
sacct -j JOB_ID --format=JobID,JobName%18,State,ReqMem,MaxRSS,AveRSS,Elapsed -P \
  | tee reproducibility_outputs/apptainer/RUN_ID/PROFILE-jhpce/slurm-resource-usage.txt
```

Replace `JOB_ID`, `RUN_ID`, and `PROFILE` with the values printed by the
runner. Keep SIF images and node-local Apptainer cache/build data outside Git.

<br><br>

## Reading build and test results

Each run writes its results to a new UTC timestamp directory, with one
subdirectory per build target (`r-portable`, `python-portable`, `r-jhpce`, or
`python-jhpce`). Use `run-info.txt` to
find the pass/fail status, then consult the log for details.

> On JHPCE, each profile directory contains the commands run, session details,
> and build/test output; the timestamp directory also contains a temporary cache clean-up log.
> The final Slurm memory and elapsed-time record is saved as `slurm-resource-usage.txt`
> in the matching `r-jhpce` or `python-jhpce` directory.

A successful run records `passed` for both the build and test status fields.
For **Docker**, `build_status=passed` means the image built; `smoke_test_status=passed`
means the checks also passed in a fresh container. For **Apptainer**, `build_status`
covers `apptainer build`, including the definition’s `%test` section.
`explicit_test_status` covers the runner’s additional `apptainer test` against
the completed SIF image. A successful **Apptainer** run also writes the SIF
image and `image.sha256`.

If a build fails, check `build.log`. If the image built but its later test
failed, check `test.log`. A missing test status usually means the run stopped
before reaching that test; for JHPCE, check the Slurm job status as well.

The test logs show the checks that ran and their output. For R, look for
the locked package count, matched versions, missing packages, and version
mismatches, followed by successful representative package imports. For Python,
look for the expected and actual Python, pySCENIC, and Arboreto versions,
`Python environment imports passed`, and `No broken requirements found.` from
`pip check`.

<br><br>

## Measured image sizes, times, and memory

These are measurements from the tested machines, not guaranteed minimum
requirements. The logs and generated images of these tests are machine-specific outputs and
remain outside Git.

### Docker on Apple Silicon (2026-11-09)

Measurements are from Docker Desktop 29.9.2 on an Apple Silicon Mac using the
`linux/amd64` emulation target. R package compilation under emulation can take
longer than native x86_64 builds.

| Image | Size | Build time | Smoke-test time |
| --- | ---: | ---: | ---: |
| R | 6.64 GB | 3,574 sec (59m 34s) | 41 sec |
| Python | 13.6 GB | 439 sec (7m 19s) | 59 sec |
| **Both** | **20.24 GB** | **About 66m 53s** | **About 1m 40s** |

Docker Desktop reported 10 CPUs and about 7.75 GiB of VM memory during these
builds; peak memory was not measured. About 35 GB free disk is a planning
estimate based on the observed image and build-cache growth, not a measured
minimum.

### JHPCE Apptainer (2026-11-09)

| Image | SIF size | Build time | Standalone test | Peak RSS | Requested memory | Total Slurm time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| R | 1.68 GiB | 5,611 sec (1h 33m 31s) | 25 sec | 4,827,528 KiB (4.61 GiB) | 10 GiB | 1h 35m 07s |
| Python | 3.80 GiB | 962 sec (16m 02s) | 7 sec | 8,506,936 KiB (8.11 GiB) | 10 GiB | 17m 33s |

In particular, the Python JHPCE build used about 8.11 GiB of its 10 GiB allocation, 
so plan additional memory on other hosts when possible.

<br><br>

## Reproducibility notes

The package lockfiles pin R and Python package versions. The base images still
use version tags instead of immutable digests, and system APT repositories are
not snapshot-pinned. These tests show that the current recipes build and pass
their smoke checks; they do not establish bit-for-bit identical images across
future rebuilds.

