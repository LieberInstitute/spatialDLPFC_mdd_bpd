#!/usr/bin/env bash
set -euo pipefail

profile=${1:?Usage: build_portable_env.sh r|python}
root=/opt/mbv

case "$profile" in
  python)
    inputs="$root/reproducibility_inputs"
    conda_root=/opt/conda
    python_env="$conda_root/envs/pyscenic"
    miniforge_version=24.7.1-2
    installer="Miniforge3-${miniforge_version}-Linux-x86_64.sh"
    url="https://github.com/conda-forge/miniforge/releases/download/${miniforge_version}/${installer}"
    curl --fail --location --silent --show-error "$url" --output /tmp/Miniforge3.sh
    bash /tmp/Miniforge3.sh -b -p "$conda_root"
    rm -f /tmp/Miniforge3.sh
    "$conda_root/bin/conda" create --yes --prefix "$python_env" \
      --file "$inputs/conda-linux-64.lock"
    (
      cd "$inputs"
      "$conda_root/bin/conda" run --no-capture-output --prefix "$python_env" \
        python -m pip install --no-deps --require-hashes \
        --requirement requirements-pip.lock
    )
    # Package caches are not part of the runtime environment. Removing them
    # keeps both Docker and Apptainer images substantially smaller.
    "$conda_root/bin/conda" clean --all --yes
    rm -rf /root/.cache/pip "$inputs"
    ;;
  r)
    target_lib="$root/renv-library"
    bootstrap_lib="$root/renv-bootstrap"
    cache_dir="$root/renv-cache"
    mkdir -p "$root/project" "$target_lib" "$bootstrap_lib" "$cache_dir"
    export R_LIBS_USER="$bootstrap_lib"
    export R_LIBS_SITE="$target_lib"
    export RENV_PATHS_CACHE="$cache_dir"
    export RENV_CONFIG_AUTO_SNAPSHOT=FALSE
    renv_url=https://cran.r-project.org/src/contrib/renv_1.3.1.tar.gz
    renv_sha256=e16281908ddf0ed6f8f6f8ad3f2cc4788af6de04c1ab223517aa1cc9acb30731
    curl --fail --location --silent --show-error "$renv_url" \
      --output /tmp/renv_1.3.1.tar.gz
    echo "$renv_sha256  /tmp/renv_1.3.1.tar.gz" | sha256sum --check --status
    Rscript --vanilla -e \
      'install.packages("/tmp/renv_1.3.1.tar.gz", repos=NULL, type="source", lib=Sys.getenv("R_LIBS_USER"))'
    rm -f /tmp/renv_1.3.1.tar.gz
    MBV_TARGET_LIB="$target_lib" Rscript --vanilla -e \
      'lib <- Sys.getenv("MBV_TARGET_LIB"); install.packages("BiocVersion", repos=c(BioCsoft="https://bioconductor.org/packages/3.20/bioc"), lib=lib); stopifnot(as.character(packageVersion("BiocVersion", lib.loc=lib)) == "3.20.0")'
    biomanager_url=https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.25.tar.gz
    biomanager_sha256=83ef8b57dcc7350fe1e8c62590a60d3ad3e1f9c8660930c946b20340a4a81e15
    curl --fail --location --silent --show-error "$biomanager_url" \
      --output /tmp/BiocManager_1.30.25.tar.gz
    echo "$biomanager_sha256  /tmp/BiocManager_1.30.25.tar.gz" | sha256sum --check --status
    MBV_TARGET_LIB="$target_lib" Rscript --vanilla -e \
      'lib <- Sys.getenv("MBV_TARGET_LIB"); install.packages("/tmp/BiocManager_1.30.25.tar.gz", repos=NULL, type="source", lib=lib); stopifnot(as.character(packageVersion("BiocManager", lib.loc=lib)) == "1.30.25")'
    rm -f /tmp/BiocManager_1.30.25.tar.gz
    cd "$root/project"
    Rscript --vanilla -e \
      'library(renv, lib.loc=Sys.getenv("R_LIBS_USER")); options(renv.install.timeout=21600); renv::restore(project=getwd(), lockfile="/opt/mbv/project/renv.lock", library="/opt/mbv/renv-library", prompt=FALSE, clean=FALSE)'
    rm -rf "$cache_dir" /root/.cache/R/renv
    ;;
  *)
    echo "Unknown profile: $profile (expected r or python)" >&2
    exit 2
    ;;
esac
