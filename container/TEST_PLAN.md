# Container validation plan

Current scope: the separate R and pySCENIC Python images support the bulk
pipeline in code directories `00` through `10`. Code directory `11_eQTL_coloc`
(tensorQTL) is deferred and may get a third image later.

| Track | Purpose | Status | Command / next action |
| --- | --- | --- | --- |
| Docker on macOS | Test user-facing Linux images via Docker Desktop, without requiring a separate Linux host | Passed on Apple Silicon with Docker Desktop 29.9.2, targeting/emulating `linux/amd64` | `bash container/scripts/run_docker_test.sh r` then `bash container/scripts/run_docker_test.sh python` |
| Portable Apptainer on clean Linux | Check general Apptainer recipes outside JHPCE | Prepared; host/access TBD, may defer | `bash container/scripts/run_portable_test.sh r` then `bash container/scripts/run_portable_test.sh python` |
| JHPCE Apptainer | Check site-specific profiles and execution on JHPCE | Waiting for explicit approval | Do not run until approved |

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
