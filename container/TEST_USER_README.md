# Test-user guide: MBv R and Python images

This guide covers **building and smoke-testing the two environment images**
for the bulk MBv pipeline in code directories `00`–`10`:

- **R image:** R 4.4.3 and the package versions in `renv.lock`.
- **Python image:** the `pyscenic_bioconda` environment, Python 3.10.14,
  pySCENIC, and Arboreto from the saved Linux x86_64 locks.

The images keep the R and Python environments separate because the code uses
them in separate jobs. Code directory `11_eQTL_coloc` (tensorQTL) is outside
this phase and may receive a third image later.

These instructions build and smoke-test environments; they do not run the MBv
analysis. The images are Docker images for Linux x86_64 (`linux/amd64`), even
when built from macOS or Windows.

## Current availability

The Docker image definitions and test runner are committed on the
`revision` branch. 

The images built in the Mac test are local to that Docker Desktop installation.
They have not been published to a container registry. The smoke-test logs and
timings are saved locally under the ignored
`reproducibility_outputs/docker/` directory.

## Build and smoke-test on Docker

Clone the project once. The lockfiles, wheel files, Dockerfiles, and test
scripts are all in that checkout; no second clone is needed.

```bash
git clone --branch revision https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd.git
cd spatialDLPFC_mdd_bpd

bash container/scripts/run_docker_test.sh r
bash container/scripts/run_docker_test.sh python
```

The runner builds and tests the R and Python images separately. It saves the
build log, smoke-test log, build duration, test duration, Docker versions,
image ID, and image size under
`reproducibility_outputs/docker/<UTC-run-id>/<profile>/`. The images are tagged
locally as `mbv-r:test` and `mbv-python:test`.

The tests check all locked R package versions and representative imports. The
Python checks verify its locked Python, pySCENIC, and Arboreto versions and
imports, then run `pip check`. They do not execute code directories `00`–`10`
end to end.

## Machine setups and compatibility

| Machine/setup | Software needed to build and smoke-test | Project clone | Compatible container form and status |
| --- | --- | --- | --- |
| Mac with Apple silicon | Docker Desktop running; Terminal/Bash | Yes; clone `revision` | Docker `linux/amd64`; **built and smoke-tested successfully** with emulation |
| Mac with Intel | [Docker Desktop for Mac](https://docs.docker.com/desktop/setup/install/mac-install/) running; Terminal/Bash | Yes; clone `revision` | Docker `linux/amd64`; target matches the machine, but this setup has not been tested |
| Windows PC with x86_64 CPU | [Docker Desktop for Windows](https://docs.docker.com/desktop/setup/install/windows-install/) using Linux containers (WSL 2 is the standard setup); run Bash commands from a WSL terminal | Yes; clone inside the WSL filesystem | Docker `linux/amd64`; target matches the machine, but Windows has not been tested |
| Linux x86_64 with Docker | Docker Engine, Git, Bash | Yes; clone `revision` | Docker `linux/amd64`; definitions target this platform, but a clean Linux Docker build has not been tested |
| Linux x86_64 with Apptainer | Apptainer, Git, Bash | Yes; clone `revision` | Portable Apptainer R/Python definitions are present; **clean Linux builds have not been tested** |
| JHPCE | JHPCE access, Slurm, and the site Apptainer module; the existing runner names `apptainer/1.5.3` | Yes; the project checkout on `revision` | Apptainer is the intended runtime, but the current JHPCE runner still targets a legacy combined image. The separate R/Python JHPCE builds are **not ready or tested**; wait for explicit approval before any JHPCE build or test |


## Image sizes, build memory, and time

The first successful build measurements are from an Apple silicon Mac using
Docker Desktop 29.9.2 and the `linux/amd64` target. This required emulation, so
the R build time is not a reliable estimate for an Intel Mac, an x86_64 PC, or
native Linux.

| Image | Final image size | Build time observed | Smoke-test time observed |
| --- | ---: | ---: | ---: |
| R | 6.64 GB | 3,574 seconds (59m 34s) | 41 seconds |
| Python | 13.6 GB | 439 seconds (7m 19s) | 59 seconds |
| **Both images** | **20.24 GB** | **About 66m 53s combined** | **About 1m 40s combined** |

Docker Desktop reported 10 CPUs and about 7.75 GiB of VM memory during the
successful Mac builds. Both images built with that allocation. Peak memory was
not measured, so treat 8 GiB assigned to Docker as the lowest configuration
demonstrated, not a measured minimum. A 16 GiB host with about 8 GiB assigned to
Docker is a practical starting point, not a tested minimum host-RAM requirement.
Allow more if available. The final image pair occupies about 20 GB. During the
observed test, host free disk space decreased by about 32 GB as images and build
layers accumulated. Plan for at least 35 GB free before building; this is a
planning estimate, and the minimum has not been measured.

The Python image is larger because its locked environment includes large
PyTorch/CUDA dependencies. It does not require a GPU to build or pass the
environment smoke checks. The first Python build failed because `pybedtools`
needed `gcc`; adding `build-essential` fixed the build. The final Python image
removes pip and Conda caches after installation.

## Portable Apptainer status

The portable definitions are:

- `container/apptainer/r-portable.def`
- `container/apptainer/python-portable.def`

On a clean Linux x86_64 host with Apptainer installed, the intended commands
are:

```bash
bash container/scripts/run_portable_test.sh r
bash container/scripts/run_portable_test.sh python
```

These builds have not been run on a clean Linux host yet. Their image sizes,
memory needs, and build times are therefore unknown. The Docker results above
do not count as portable Apptainer test results.
