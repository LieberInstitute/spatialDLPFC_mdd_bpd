# MBv pipeline containers

For this container phase, the target is the bulk of the pipeline in code
directories `00` through `10`, which uses R and Python in separate jobs. This
repository therefore builds separate R and pySCENIC Python images; the R image
does not put the pySCENIC Conda interpreter on `PATH`. A user needs one clone
of the MBv repository.
The lockfiles, wheels, build definitions, and verification scripts are included
in that checkout. The environment smoke tests do not require the study data.
The tensorQTL work in code directory `11_eQTL_coloc` is outside this phase and
may receive its own image later.

## Docker on macOS or Windows

Docker Desktop runs Linux containers inside its Linux VM. It lets Mac and
Windows users run the Linux environments without managing a separate Linux
machine. The current Conda lock is `linux-64`, so both images are explicitly
built as `linux/amd64`. On Apple Silicon this uses architecture emulation, and
the R package compilation may be slow. Docker documents that emulated builds,
especially compilation, can take considerably longer than native builds.

With Docker Desktop started, clone this repository and build/test each image:

```bash
git clone --branch revision https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd.git
cd spatialDLPFC_mdd_bpd
bash container/scripts/run_docker_test.sh r
bash container/scripts/run_docker_test.sh python
```

The runner builds each image, checks locked package versions and representative
imports, then runs the checks again in a fresh container. It records the exact
command, platform, Git commit, Docker versions, build duration, smoke-test
duration, image ID and size, and logs under
`reproducibility_outputs/docker/<UTC-run-id>/<profile>/`. Outputs are ignored
by Git so image files and machine-specific logs stay out of repository history.
The images remain local to that Docker installation; publishing them to a
registry is a separate step.

The included Dockerfiles currently target Linux x86_64 to match the saved
Python lock. Supporting native `linux/arm64` needs a separate Conda lock and
validation. Docker Desktop on Apple Silicon supports emulation for non-native
platform builds; see [Docker's multi-platform build documentation](https://docs.docker.com/build/building/multi-platform/).

### Initial Apple Silicon build measurements

On an Apple Silicon Mac using Docker Desktop 29.9.2 and the `linux/amd64`
emulation target, the first successful R image build took 3,574 seconds
(59 minutes 34 seconds), followed by a 41-second smoke test. The Python image
build took 439 seconds (7 minutes 19 seconds), followed by a 59-second smoke
test. The resulting images were 6.64 GB (R) and 13.6 GB (Python). The Python
image includes the locked CUDA/PyTorch packages even when the host has no
NVIDIA GPU. Allow at least about 20 GB for the final pair of images, plus
temporary build cache. Treat these as measurements for this Mac and these
package sources, not estimates for Intel Macs, PCs, or native Linux.

The first Python attempt failed because `pybedtools` needed to compile but the
image lacked `gcc`; adding `build-essential` fixed it. That compiler is included
in both portable Python profiles. The final Python image removes pip and Conda
package caches after installation, then verifies imports and package consistency.

## Portable Apptainer on Linux

The portable Apptainer definitions are `apptainer/r-portable.def` and
`apptainer/python-portable.def`. On a clean Linux x86_64 machine with Apptainer
installed, clone the repository and run:

```bash
git clone --branch revision https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd.git
cd spatialDLPFC_mdd_bpd
bash container/scripts/run_portable_test.sh r
bash container/scripts/run_portable_test.sh python
```

The runner defaults to `apptainer build --fakeroot`, which requires fakeroot
support configured in Apptainer. If the host has root access but no fakeroot
support, run it as root with `APPTAINER_BUILD_MODE=root`. Logs, build and test
durations, image checksums, and metadata are saved under
`reproducibility_outputs/apptainer/<UTC-run-id>/<profile>/`.

## JHPCE Apptainer

The JHPCE workflow remains a distinct site-specific track. Do not run
`container/scripts/run_jhpce_apptainer_test.sh` until that phase is approved.
Its existing combined definition is retained for that workflow while the
portable profiles are validated.

## Reproducibility limits

The smoke tests compare all R package versions to `renv.lock`, exercise
representative R imports, check Python/pySCENIC/Arboreto versions and imports,
and run `pip check`. They do not yet run the code-directory `00` through `10`
analysis end to end or validate external study data.

The base images use version tags rather than immutable digests, and Debian APT
repositories are not yet snapshot-pinned. The first tests will establish
whether the dependency recipes work on each target; pinning base digests and
system-package snapshots remains necessary before claiming bit-for-bit image
rebuilds.
