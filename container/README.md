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

## Build and distribution directory layout

| Order | Associated directory in `container/` | Description | When to run |
| --- | --- | --- | --- |
| Stage 1 | `01_environment/` | Holds the frozen Python environment inputs. The R lockfile and `renv/` activation directory remain at the project root; see the [environment README](01_environment/README.md) for updating the R lockfile. | Only necessary if future R/Python packages are installed requiring new image builds. |
| Stage 2 | `02_image_build_test/` | Holds Docker and Apptainer definitions, build and smoke-test scripts, and the build/test plan; see [build/test README](02_image_build_test_README.md) for instructions. | When building images or checking that their environments install and pass smoke tests. |
| Stage 3 | `03_image_code_check/` | Holds representative project code checks and instructions for staging their inputs; see [code check README](03_image_code_check/README.md) for instructions. | After building an image, to check it with representative project code. |
| Stage 4 | `04_image_distribution/` | Holds the [distribution plan](04_image_distribution/DISTRIBUTION_PLAN.md). | When preparing to publish or share images. |
| — | `outputs/` | Holds generated build and test artifacts and logs separately from source files. | Not run, contains the outputs for Stage 2-3. |

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
| R | `container/02_image_build_test/apptainer/r-portable.def` | `container/02_image_build_test/apptainer/r-jhpce.def` |
| Python | `container/02_image_build_test/apptainer/python-portable.def` | `container/02_image_build_test/apptainer/python-jhpce.def` |

The portable definitions (**yet to be tested**) use normal Linux build behavior.
The JHPCE definitions apply site-specific workarounds required by its 
rootless/fakeroot build environment: APT runs with its sandbox user set to root, 
archive extraction uses `TAR_OPTIONS=--no-same-owner`, and staged verifier and R
lockfile permissions are normalized for runtime checks. Both profiles call
shared build scripts, so package installation logic stays aligned while the
JHPCE adjustments remain explicit in the JHPCE profile.

