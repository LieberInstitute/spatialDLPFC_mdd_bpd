
# Stage 2: Build and test MBv container images

This guide explains how to build and smoke-test the separate R and Python
images for the bulk MBv pipeline environment in code directories `00`–`10`.
It covers container setup only, not implementation or execution of the
analysis workflow.

## Tested environments and availability

The JHPCE Apptainer and MacOS Docker R and Python images were built and smoke-tested 
successfully on 2026-11-12 (see details [below](#measured-image-sizes-times-and-memory)). 
R builds consistently take over an hour, so plan accordingly.

Docker images built with the commands below are local to that Docker
installation; they have not been published to a container registry.

<br>

Use this table to choose the build instructions best-suited for your setup.

| Machine/setup | Software needed to build and smoke-test | Project clone | Compatible image and current test status |
| --- | --- | --- | --- |
| Mac with Apple Silicon | Docker Desktop running; Terminal/Bash | Yes; clone `revision` | Docker `linux/amd64`; latest R and Python builds and checks passed on 2026-11-12 |
| Mac with Intel | Docker Desktop for Mac; Terminal/Bash | Yes; clone `revision` | Docker `linux/amd64`; matches CPU architecture, but this setup has not been tested |
| Windows PC with x86_64 CPU | Docker Desktop using Linux containers; WSL 2 and Bash terminal recommended | Yes; clone inside WSL filesystem | Docker `linux/amd64`; matches CPU architecture, but this setup has not been tested |
| Linux x86_64 with Docker | Docker Engine, Git, Bash | Yes; clone `revision` | Docker `linux/amd64`; clean Linux Docker build has not been tested |
| Linux x86_64 with Apptainer | Apptainer, Git, Bash; configured fakeroot or root access | Yes; clone `revision` | Portable Apptainer definitions are prepared; clean Linux builds have not been tested |
| JHPCE | JHPCE access, Slurm, and Apptainer module `apptainer/1.5.3` | Yes; checkout on `revision` | Separate named JHPCE Apptainer definitions; latest R and Python builds and checks passed on 2026-11-12 |

The Python lock targets Linux x86_64, so all Docker builds use
`linux/amd64`. Docker Desktop on Apple Silicon emulates that architecture;
native `linux/arm64` images would need a separate Conda lock and validation.

<br>

## Build and test steps

Briefly, you need to clone the project and then use bash to build the container image(s)
and test the builds. These tests verify the environment, not the full pipeline.
The build logs and test results are stored in `container/outputs/image_build_test`
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
# navigate to cloned directory
cd spatialDLPFC_mdd_bpd
# run build/test
bash container/02_image_build_test/scripts/run_docker_build_smoke_test.sh r
bash container/02_image_build_test/scripts/run_docker_build_smoke_test.sh python
```

The runner builds each image, checks its locked versions and representative
imports, and repeats the checks in a fresh container. Docker images are tagged
locally as `mbv-r:test` and `mbv-python:test`. The runner saves commands, logs,
build and test durations, Docker versions, image ID, and image size under
`container/outputs/image_build_test/docker/<UTC-run-id>/<profile>/`.

On macOS, Docker Desktop must be running and have enough VM memory and disk.
On Windows, use Linux containers and run the commands from WSL 2.

<br>

### Build Option 2: Apptainer images on Linux

***These portable definitions have not yet been tested on a clean Linux host.***

Use a Linux x86_64 host with Apptainer installed, then run R and Python as
separate builds:

```bash
# navigate to cloned directory
cd spatialDLPFC_mdd_bpd
# run build/test
bash container/02_image_build_test/scripts/run_apptainer_portable_build_smoke_test.sh r
bash container/02_image_build_test/scripts/run_apptainer_portable_build_smoke_test.sh python
```

The runner defaults to `apptainer build --fakeroot`, which requires fakeroot
support configured in Apptainer. If the host has root access but no fakeroot
support, run the commands as root with `APPTAINER_BUILD_MODE=root`. The
portable definitions use normal Linux APT behavior and do not include the
JHPCE-specific settings. Logs, build and test durations, checksums, and
metadata go under `container/outputs/image_build_test/apptainer/<UTC-run-id>/<profile>-portable/`.

<br>

### Build Option 3: Apptainer images on JHPCE

The JHPCE definitions are `container/02_image_build_test/apptainer/r-jhpce.def` and
`container/02_image_build_test/apptainer/python-jhpce.def`. Use the login node to request a
separate 10 GiB compute allocation for each image. The runner loads the
Apptainer module, builds the selected image, and performs its smoke test.

```bash
# navigate to cloned directory
cd spatialDLPFC_mdd_bpd

# make dir for batch job outputs
mkdir -p container/outputs/image_build_test/apptainer/batch-logs

# Option: Run build/test as batch job for both images
sbatch \
  --job-name=mbv-all-jhpce \
  --mem=10G \
  --output="container/outputs/image_build_test/apptainer/batch-logs/all-%j.out" \
  --wrap='bash container/02_image_build_test/scripts/run_apptainer_jhpce_build_smoke_test.sh all'

# Option: Run build/test as batch job for single profile (python example)
sbatch \
  --job-name=mbv-python-jhpce \
  --mem=10G \
  --output="container/outputs/image_build_test/apptainer/batch-logs/python-%j.out" \
  --wrap='bash container/02_image_build_test/scripts/run_apptainer_jhpce_build_smoke_test.sh python'

# Option: run build/test as interactive job (single profile, R example)
srun --mem 10G --pty bash
export MBV_JHPCE_MEM_REQUEST=10G
bash container/02_image_build_test/scripts/run_apptainer_jhpce_build_smoke_test.sh r
## job runs
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
`container/outputs/image_build_test/apptainer/<UTC-run-id>/<profile>-jhpce/`. Use the
`session-info.txt` file to batch the UTC run ID to the job ID.

<br><br>

## Reading build and test results

Each run writes its results to a new UTC timestamp directory, with one
subdirectory per build target (`r-portable`, `python-portable`, `r-jhpce`, or
`python-jhpce`). Use `run-info.txt` to
find the pass/fail status, then consult the log for details.

> On JHPCE, each profile directory contains the commands run, session details,
> and build/test output; the timestamp directory also contains a temporary cache clean-up log.
> If using a batch job, the log saved in `container/outputs/image_build_test/apptainer/batch-logs`
> for the specific job ID will also report the test status and results path.

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

Apptainer `.sif` images are stored in `container/outputs/image_build_test/apptainer/` profile- and run-specific subdirectories. These files are not tracked by Git due to their large size.

<br><br>

## Measured image sizes, times, and memory

These are measurements from the tested machines, not guaranteed minimum
requirements. The logs and generated images of these tests are machine-specific outputs and
remain outside Git.

### Docker on Apple Silicon (measurements through 2026-11-12)

Measurements are from Docker Desktop 29.9.2 on an Apple Silicon Mac using the
`linux/amd64` emulation target. R package compilation under emulation can take
longer than native x86_64 builds. Docker's image size reports the image's 
uncompressed contents.

| Image | Size | Build time | Smoke-test time |
| --- | ---: | ---: | ---: |
| R | 6.33 GiB (6.80 GB) | 3,332 sec (55m 32s) | 40 sec |
| Python | 12.67 GiB (13.60 GB) | 425 sec (7m 05s) | 50 sec |
| **Both** | **19.00 GiB (20.40 GB)** | **About 1hr 3min** | **About 1.5 min** |

About 35 GB free disk is a planning estimate based on the observed
image and build-cache growth, not a measured minimum.

### JHPCE Apptainer (2026-11-12)

Apptainer's image size is the `.sif` file's on-disk size which is compressed.

| Image | SIF size | Build time | Standalone test | Peak RSS | Requested memory |
| --- | ---: | ---: | ---: | ---: | ---: |
| R | 1.73 GiB (1.86 GB) | 6,428 sec (1h 47m 08s) | 28 sec | Not recorded | 10G |
| Python | 3.79 GiB (4.07 GB) | 977 sec (16m 17s) | 7 sec | Not recorded | 10G |

Build times
may vary depending on node assignment.

<br><br>

## Reproducibility notes

The package lockfiles pin R and Python package versions. The base images still
use version tags instead of immutable digests, and system APT repositories are
not snapshot-pinned. These tests show that the current recipes build and pass
their smoke checks; they do not establish bit-for-bit identical images across
future rebuilds.

