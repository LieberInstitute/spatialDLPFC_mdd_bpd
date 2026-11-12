# Image build and test scripts

This directory contains scripts used by Docker and Apptainer image builds and
their environment smoke tests. Project-level checks that run AUCell or construct
a `SpatialExperiment` are in `container/03_image_code_check/`.

| File | Role | Called from |
| --- | --- | --- |
| `install_locked_environment.sh` | Installs the locked R or Python environment into the image. This is the shared package installation helper. | Dockerfiles and `prepare_apptainer_environment.sh` |
| `prepare_apptainer_environment.sh` | Installs system packages, applies the portable or JHPCE Apptainer build settings, then calls the shared environment installer. | Apptainer definition `%post` sections |
| `verify_r_environment.R` | Compares installed R package versions with `renv.lock`. | R Docker build and R Apptainer `%test` |
| `verify_python_environment.py` | Checks Python, pySCENIC, and Arboreto versions and representative imports. | Python Docker build and Python Apptainer `%test` |
| `run_docker_build_smoke_test.sh` | Builds one Docker profile and runs its environment smoke tests on the host. | User command from the project root |
| `run_apptainer_portable_build_smoke_test.sh` | Builds one portable Apptainer profile and runs its environment smoke tests on the host. | User command on Linux |
| `run_apptainer_jhpce_build_smoke_test.sh` | Builds one or both JHPCE Apptainer profiles and runs their environment smoke tests in a Slurm allocation. | User command inside a JHPCE compute allocation |

The Docker and Apptainer definitions select which profile is built. The
verifier scripts are copied into the images under `/opt/mbv/scripts/`; the
`run_*` scripts remain host-side orchestration tools.
