# Container validation plan

Current scope: the separate R and pySCENIC Python images support the bulk
pipeline in code directories `00` through `10`. Code directory `11_eQTL_coloc`
(tensorQTL) is deferred and may get a third image later.

| Track | Purpose | Status | Command / next action |
| --- | --- | --- | --- |
| Docker on macOS | Test user-facing Linux images via Docker Desktop, without requiring a separate Linux host | Passed on Apple Silicon with Docker Desktop 29.9.2, targeting/emulating `linux/amd64` | `bash container/scripts/run_docker_test.sh r` then `bash container/scripts/run_docker_test.sh python` |
| Portable Apptainer on clean Linux | Check general Apptainer recipes outside JHPCE | Prepared; host/access TBD, may defer | `bash container/scripts/run_portable_test.sh r` then `bash container/scripts/run_portable_test.sh python` |
| JHPCE Apptainer | Check separate site-specific R and Python builds on JHPCE | Passed on 2026-11-09; both images built and passed separate runtime checks in 10 GB compute allocations | Complete; see measured results below |

The first two tracks use `linux/amd64` because the checked-in Conda explicit
lock targets Linux x86_64. Docker Desktop on Apple Silicon emulates that target;
record Mac build durations separately from native Linux timings. Keep all
images and logs under the ignored `reproducibility_outputs/` tree.


## Apple Silicon Docker results (2026-11-09)

| Image | Build | Fresh-container smoke test | Image size | Result directory |
| --- | ---: | ---: | ---: | --- |
| R | 3,574 sec (59m 34s) | 41 sec | 6.64 GB | `reproducibility_outputs/docker/20261109T154726Z/r/` |
| Python | 439 sec (7m 19s) | 59 sec | 13.6 GB | `reproducibility_outputs/docker/20261109T170322Z/python/` |

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

## JHPCE Apptainer results (2026-11-09)

| Image | Build time | Standalone test | SIF size | Peak RSS (10 GB requested) | Result directory |
| --- | ---: | ---: | ---: | ---: | --- |
| R | 5,611 sec (1h 33m 31s) | 25 sec | 1.68 GiB | 4,827,528 KiB (4.61 GiB) | `reproducibility_outputs/apptainer/20261109T214031Z/jhpce/r/` |
| Python | 962 sec (16m 02s) | 7 sec | 3.80 GiB | 8,506,936 KiB (8.11 GiB) | `reproducibility_outputs/apptainer/20261109T231706Z/jhpce/python/` |

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
`reproducibility_outputs/apptainer/20261109T191626Z/jhpce/r/`.
