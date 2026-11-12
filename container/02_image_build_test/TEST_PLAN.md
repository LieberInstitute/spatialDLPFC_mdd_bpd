# Container validation plan

Current scope: the separate R and pySCENIC Python images support the bulk
pipeline in code directories `00` through `10`. Code directory `11_eQTL_coloc`
(tensorQTL) is deferred and may get a third image later.

| Track | Purpose | Status | Command / next action |
| --- | --- | --- | --- |
| Docker on macOS | Test user-facing Linux images via Docker Desktop, without requiring a separate Linux host | Passed on Apple Silicon with Docker Desktop 29.9.2, targeting/emulating `linux/amd64` | `bash container/02_image_build_test/scripts/run_docker_build_smoke_test.sh r` then `bash container/02_image_build_test/scripts/run_docker_build_smoke_test.sh python` |
| Portable Apptainer on clean Linux | Check general Apptainer recipes outside JHPCE | Prepared; host/access TBD, may defer | `bash container/02_image_build_test/scripts/run_apptainer_portable_build_smoke_test.sh r` then `bash container/02_image_build_test/scripts/run_apptainer_portable_build_smoke_test.sh python` |
| JHPCE Apptainer | Check separate site-specific R and Python builds on JHPCE | Latest R and Python builds and explicit tests passed on 2026-11-12 in one 10G Slurm batch job | Complete; see measured results below |

The first two tracks use `linux/amd64` because the checked-in Conda explicit
lock targets Linux x86_64. Docker Desktop on Apple Silicon emulates that target;
record Mac build durations separately from native Linux timings. Keep generated
images and logs under `container/outputs/image_build_test/`; `.sif` images are
ignored by Git, while logs and metadata may be tracked.
New Apptainer runs use `container/outputs/image_build_test/apptainer/<UTC-run-id>/<profile>-<site>/`,
for example `r-portable` or `python-jhpce`. The completed JHPCE result paths
below retain their original layout because those records were created before
the output convention changed.


## Apple Silicon Docker results (2026-11-09)

| Image | Build | Fresh-container smoke test | Image size | Result directory |
| --- | ---: | ---: | ---: | --- |
| R | 3,574 sec (59m 34s) | 41 sec | 6.64 GB | `container/outputs/image_build_test/docker/20261109T154726Z/r/` |
| Python | 439 sec (7m 19s) | 59 sec | 13.6 GB | `container/outputs/image_build_test/docker/20261109T170322Z/python/` |

R verified all 415 locked package versions with no missing packages or
mismatches, then loaded representative spatial and single-cell libraries.
Python verified Python 3.10.14, pySCENIC 0.12.1+11.g06bafba, Arboreto 0.1.6,
their imports, and `pip check`.

The first Python attempt failed after 189 sec because `pybedtools` needed
`gcc`. Adding `build-essential` to both portable Python definitions fixed the
build. The final Python image also removes pip and Conda caches; its size
decreased from 20.9 GB to 13.6 GB and the smoke test still passed. An initial R
attempt stopped before downloading the base image because Docker Desktop's
credential helper was not on the calling shell's `PATH`; the runner now adds
Docker Desktop's resource directory on macOS.

## Apple Silicon Docker Python rebuild (2026-11-12)

The restructured Python Docker build and smoke test passed on Docker Desktop
29.9.2, targeting `linux/amd64` under Apple Silicon emulation. The image size
was 13,604,890,545 bytes (12.67 GiB / 13.60 GB). Build time was 425 seconds
(7m 05s), and the fresh-container smoke test took 50 seconds. The run record,
build log, and smoke-test log are in
`container/outputs/image_build_test/docker/20261112T164344Z/python/`.
The image was built from commit `c54dd52a` and has Docker image ID
`sha256:94c34217a023d8f94de4205363d1cc18b8797f3d8e863789898e10768c92005e`.

The Docker AUCell code check also completed on the rebuilt Python image using
the first 10% (2,126 of 21,265 cells). Its log and output loom are under
`container/outputs/image_code_check/docker/`.

## Apple Silicon Docker R rebuild (2026-11-12)

The restructured R Docker image built and passed its smoke test on Docker
Desktop 29.9.2, targeting `linux/amd64` under Apple Silicon emulation. The image
size was 6,798,736,903 bytes (6.33 GiB / 6.80 GB). Build time was 3,332 seconds
(55m 32s), and the fresh-container smoke test took 40 seconds. The run record,
build log, and smoke-test log are in
`container/outputs/image_build_test/docker/20261112T165849Z/r/`.
The image was built from commit `83c4597b` and has Docker image ID
`sha256:8b1779430c71e450134a20562ec10aea2a27a577e1beb8fa5d2a985b2014360d`.

The R verifier matched all 425 packages in `renv.lock` with no missing packages
or version mismatches. Representative spatial and single-cell packages loaded
successfully. During the build, Docker Desktop was configured with 10 CPUs and
8 GiB of VM memory; peak memory use was not recorded.

The Docker R code check downloaded the 4.2 GB `GSE333974_RAW.tar` archive,
created an SPE from three samples, reached `*** TEST END ***`, and wrote
`spe_n3_se.rds` and `spe_n3_assays.h5`. Its log and outputs are under
`container/outputs/image_code_check/docker/`.

## JHPCE Apptainer results (2026-11-09)

| Image | Build time | Standalone test | SIF size | Peak RSS (10 GB requested) | Result directory |
| --- | ---: | ---: | ---: | ---: | --- |
| R | 5,611 sec (1h 33m 31s) | 25 sec | 1.68 GiB | 4,827,528 KiB (4.61 GiB) | `container/outputs/image_build_test/apptainer/20261109T214031Z/jhpce/r/` |
| Python | 962 sec (16m 02s) | 7 sec | 3.80 GiB | 8,506,936 KiB (8.11 GiB) | `container/outputs/image_build_test/apptainer/20261109T231706Z/jhpce/python/` |

R matched all 415 package versions in `renv.lock`, with no missing packages or
version mismatches, and loaded the representative spatial and single-cell
libraries. Python verified Python 3.10.14, pySCENIC 0.12.1+11.g06bafba, and
Arboreto 0.1.6; required imports passed and `pip check` reported no broken
requirements. Total Slurm job times were 1h 35m 07s for R and 17m 33s for
Python. Both runs used separate `srun --mem 10G --pty bash` allocations.

The first R attempt built its SIF but the standalone test could not read
`/opt/mbv/project/renv.lock`. The shared build helper now sets readable
permissions on the project directory and copied lockfile; the subsequent R
build and test passed. The failed run remains preserved in
`container/outputs/image_build_test/apptainer/20261109T191626Z/jhpce/r/`.

## JHPCE Apptainer rebuild (2026-11-12)

Both restructured JHPCE images built and passed `apptainer test` in batch job
`37038242` on compute node `compute-168.cm.cluster`, using Apptainer
`1.5.3+342-g79cf6087c`. The job requested `10G` of memory. Peak RSS was not
recorded in the saved run files.

| Image | Build time | Explicit test | SIF size | SHA-256 |
| --- | ---: | ---: | ---: | --- |
| R | 6,428 sec (1h 47m 08s) | 28 sec | 1,858,621,440 bytes (1.73 GiB / 1.86 GB) | `7242e59319c90b438c6f06305c9f9086ebbb0b628d821cf5597a88540fdecfd2` |
| Python | 977 sec (16m 17s) | 7 sec | 4,074,225,664 bytes (3.79 GiB / 4.07 GB) | `49003df1e082aa54ab072897f95e2a3964548d4961e862eec419d228d1017e6b` |

R matched all 425 package versions in `renv.lock`, with no missing packages or
version mismatches, and loaded the representative spatial and single-cell
libraries. Python verified Python 3.10.14, pySCENIC 0.12.1+11.g06bafba, and
Arboreto 0.1.6; imports passed and `pip check` reported no broken requirements.
The run records, build/test logs, definitions used, and checksums are under
`container/outputs/image_build_test/apptainer/20261112T171336Z/` on JHPCE.
