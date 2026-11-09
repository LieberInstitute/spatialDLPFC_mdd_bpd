FROM rocker/r-ver:4.4.3

LABEL org.opencontainers.image.title="MBv R environment" \
      org.opencontainers.image.version="0.1-test" \
      org.opencontainers.image.description="R renv environment for the MBv pipeline"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran cmake pkg-config git curl ca-certificates \
        libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libssh2-1-dev \
        libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libcairo2-dev libharfbuzz-dev libfribidi-dev libxt-dev libx11-dev \
        gdal-bin libgdal-dev libgeos-dev libproj-dev libudunits2-dev libhdf5-dev \
        libnetcdf-dev libglpk-dev libgsl-dev libfftw3-dev libopenblas-dev \
        liblapack-dev libblas-dev libpq-dev libmagick++-dev libzmq3-dev \
        libicu-dev libarchive-dev libtbb-dev libabsl-dev pandoc \
        zlib1g-dev libbz2-dev liblzma-dev libpcre2-dev \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/mbv/project /opt/mbv/scripts

COPY renv.lock /opt/mbv/project/renv.lock
COPY container/scripts/verify_r_lockfile.R /opt/mbv/scripts/verify_r_lockfile.R
COPY container/scripts/build_portable_env.sh /opt/mbv/scripts/build_portable_env.sh

ENV MBV_R_LIBRARY=/opt/mbv/renv-library \
    R_LIBS_SITE=/opt/mbv/renv-library \
    R_LIBS_USER=/opt/mbv/renv-bootstrap

RUN chmod a+r /opt/mbv/scripts/verify_r_lockfile.R \
    && bash /opt/mbv/scripts/build_portable_env.sh r \
    && Rscript --vanilla /opt/mbv/scripts/verify_r_lockfile.R \
        /opt/mbv/renv-library /opt/mbv/project/renv.lock \
    && Rscript --vanilla -e 'library(sf); library(s2); library(Seurat); library(SingleCellExperiment); library(spatialLIBD)'

WORKDIR /work
CMD ["Rscript", "--version"]
