#!/usr/bin/env bash
set -euo pipefail
# Package source archives record upstream owners that are not mapped inside
# JHPCE's rootless Apptainer user namespace. Ownership is irrelevant for the
# installed files, so prevent GNU tar from trying to restore it during R
# package unpacking.
export TAR_OPTIONS=--no-same-owner

ROOT=/opt/mbv
INPUTS="$ROOT/reproducibility_inputs"
TARGET_LIB="$ROOT/renv-library"
BOOTSTRAP_LIB="$ROOT/renv-bootstrap"
CACHE_DIR="$ROOT/renv-cache"
CONDA_ROOT=/opt/conda
PYTHON_ENV="$CONDA_ROOT/envs/pyscenic"
MINIFORGE_VERSION=24.7.1-2
MINIFORGE_INSTALLER="Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh"
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_INSTALLER}"

mkdir -p "$ROOT/project" "$TARGET_LIB" "$BOOTSTRAP_LIB" "$CACHE_DIR"

# Install a versioned Conda bootstrapper, then use the saved explicit lock and
# hash-checked pip requirements to create the pySCENIC environment.
curl --fail --location --silent --show-error "$MINIFORGE_URL" \
    --output /tmp/Miniforge3.sh
bash /tmp/Miniforge3.sh -b -p "$CONDA_ROOT"
rm -f /tmp/Miniforge3.sh
"$CONDA_ROOT/bin/conda" create --yes --prefix "$PYTHON_ENV" \
    --file "$INPUTS/conda-linux-64.lock"
(
    cd "$INPUTS"
    "$CONDA_ROOT/bin/conda" run --no-capture-output --prefix "$PYTHON_ENV" \
        python -m pip install --no-deps --require-hashes \
        --requirement requirements-pip.lock
)

# Restore the R lockfile in an isolated library. BiocVersion is seeded from
# Bioconductor 3.20 first because the known-working JHPCE restore required it.
export R_LIBS_USER="$BOOTSTRAP_LIB"
export R_LIBS_SITE="$TARGET_LIB"
export RENV_PATHS_CACHE="$CACHE_DIR"
export RENV_CONFIG_AUTO_SNAPSHOT=FALSE
RENV_BOOTSTRAP_URL="https://cran.r-project.org/src/contrib/renv_1.3.1.tar.gz"
RENV_BOOTSTRAP_SHA256=e16281908ddf0ed6f8f6f8ad3f2cc4788af6de04c1ab223517aa1cc9acb30731
curl --fail --location --silent --show-error "$RENV_BOOTSTRAP_URL" \
    --output /tmp/renv_1.3.1.tar.gz
echo "$RENV_BOOTSTRAP_SHA256  /tmp/renv_1.3.1.tar.gz" | sha256sum --check --status
Rscript --vanilla -e \
    'install.packages("/tmp/renv_1.3.1.tar.gz", repos=NULL, type="source", lib=Sys.getenv("R_LIBS_USER"))'
rm -f /tmp/renv_1.3.1.tar.gz
MBV_TARGET_LIB="$TARGET_LIB" Rscript --vanilla -e \
    'lib <- Sys.getenv("MBV_TARGET_LIB"); install.packages("BiocVersion", repos=c(BioCsoft="https://bioconductor.org/packages/3.20/bioc"), lib=lib); stopifnot(as.character(packageVersion("BiocVersion", lib.loc=lib)) == "3.20.0")'
BIOMANAGER_URL="https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.25.tar.gz"
BIOMANAGER_SHA256=83ef8b57dcc7350fe1e8c62590a60d3ad3e1f9c8660930c946b20340a4a81e15
curl --fail --location --silent --show-error "$BIOMANAGER_URL" \
    --output /tmp/BiocManager_1.30.25.tar.gz
echo "$BIOMANAGER_SHA256  /tmp/BiocManager_1.30.25.tar.gz" | sha256sum --check --status
MBV_TARGET_LIB="$TARGET_LIB" Rscript --vanilla -e \
    'lib <- Sys.getenv("MBV_TARGET_LIB"); install.packages("/tmp/BiocManager_1.30.25.tar.gz", repos=NULL, type="source", lib=lib); stopifnot(as.character(packageVersion("BiocManager", lib.loc=lib)) == "1.30.25")'
rm -f /tmp/BiocManager_1.30.25.tar.gz
cd "$ROOT/project"
Rscript --vanilla -e \
    'library(renv, lib.loc=Sys.getenv("R_LIBS_USER")); options(renv.install.timeout=21600); renv::restore(project=getwd(), lockfile="renv.lock", library="/opt/mbv/renv-library", prompt=FALSE, clean=FALSE)'

# The installed library is retained; the temporary renv download cache is not.
rm -rf "$CACHE_DIR" /root/.cache/R/renv
