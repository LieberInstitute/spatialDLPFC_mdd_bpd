#!/usr/bin/env bash
set -euo pipefail

profile=${1:?Usage: prepare_apptainer_environment.sh r|python portable|jhpce}
target=${2:?Usage: prepare_apptainer_environment.sh r|python portable|jhpce}

case "$profile" in
  r)
    apt_packages=(
      build-essential gfortran cmake pkg-config git curl ca-certificates
      libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libssh2-1-dev
      libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
      libcairo2-dev libharfbuzz-dev libfribidi-dev libxt-dev libx11-dev
      gdal-bin libgdal-dev libgeos-dev libproj-dev libudunits2-dev libhdf5-dev
      libnetcdf-dev libglpk-dev libgsl-dev libfftw3-dev libopenblas-dev
      liblapack-dev libblas-dev libpq-dev libmagick++-dev libzmq3-dev
      libicu-dev libarchive-dev libtbb-dev libabsl-dev pandoc
      zlib1g-dev libbz2-dev liblzma-dev libpcre2-dev
    )
    verifier=/opt/mbv/scripts/verify_r_environment.R
    ;;
  python)
    apt_packages=(build-essential ca-certificates curl bzip2 procps)
    verifier=/opt/mbv/scripts/verify_python_environment.py
    ;;
  *) echo "Profile must be r or python" >&2; exit 2 ;;
esac

case "$target" in
  portable)
    apt_command=(apt-get)
    unset TAR_OPTIONS || true
    ;;
  jhpce)
    # JHPCE's fakeroot user namespace cannot use APT's default _apt sandbox.
    apt_command=(apt-get -o APT::Sandbox::User=root)
    # Source archives may contain owners unavailable inside JHPCE fakeroot.
    export TAR_OPTIONS=--no-same-owner
    ;;
  *) echo "Target must be portable or jhpce" >&2; exit 2 ;;
esac

export DEBIAN_FRONTEND=noninteractive
"${apt_command[@]}" update
"${apt_command[@]}" install -y --no-install-recommends "${apt_packages[@]}"
rm -rf /var/lib/apt/lists/*

# Build staging can assign restrictive host group permissions to source files.
chmod 644 "$verifier"
if [[ $profile == r ]]; then
  chmod 755 /opt/mbv/project
  chmod 644 /opt/mbv/project/renv.lock
fi
bash /opt/mbv/scripts/install_locked_environment.sh "$profile"
