# syntax=docker/dockerfile:1.7

FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="MBv Python environment" \
      org.opencontainers.image.version="0.1-test" \
      org.opencontainers.image.description="pySCENIC Conda environment for the MBv pipeline"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential bash ca-certificates curl bzip2 procps \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/mbv/scripts

COPY container/02_image_build_test/scripts/verify_python_environment.py /opt/mbv/scripts/verify_python_environment.py
COPY container/02_image_build_test/scripts/install_locked_environment.sh /opt/mbv/scripts/install_locked_environment.sh

ENV PATH=/opt/conda/envs/pyscenic/bin:/opt/conda/bin:$PATH

# Mount lockfiles and wheels only during installation so they don't add an image layer.
RUN --mount=type=bind,source=container/01_environment/python,target=/tmp/mbv-build-inputs/python,readonly \
    export MBV_BUILD_INPUTS=/tmp/mbv-build-inputs/python \
           MBV_BUILD_INPUTS_EPHEMERAL=1 \
    && chmod a+r /opt/mbv/scripts/verify_python_environment.py \
    && bash /opt/mbv/scripts/install_locked_environment.sh python \
    && python /opt/mbv/scripts/verify_python_environment.py \
    && python -m pip check

# Confirm the temporary build inputs were not persisted in the image.
RUN test ! -e /opt/mbv/build_inputs/python/conda-linux-64.lock \
    && test ! -e /opt/mbv/build_inputs/python/requirements-pip.lock \
    && test ! -e /opt/mbv/build_inputs/python/wheels \
    && test ! -e /tmp/mbv-build-inputs/python/conda-linux-64.lock \
    && test ! -e /tmp/mbv-build-inputs/python/requirements-pip.lock \
    && test ! -e /tmp/mbv-build-inputs/python/wheels

WORKDIR /work
CMD ["python", "--version"]
