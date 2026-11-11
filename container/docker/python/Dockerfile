# syntax=docker/dockerfile:1.7

FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="MBv Python environment" \
      org.opencontainers.image.version="0.1-test" \
      org.opencontainers.image.description="pySCENIC Conda environment for the MBv pipeline"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential bash ca-certificates curl bzip2 procps \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/mbv/scripts

COPY container/scripts/verify_python_env.py /opt/mbv/scripts/verify_python_env.py
COPY container/scripts/build_portable_env.sh /opt/mbv/scripts/build_portable_env.sh

ENV PATH=/opt/conda/envs/pyscenic/bin:/opt/conda/bin:$PATH

# Mount lockfiles and wheels only during installation so they don't add an image layer.
RUN --mount=type=bind,source=reproducibility_inputs,target=/tmp/mbv-reproducibility_inputs,readonly \
    export MBV_REPRODUCIBILITY_INPUTS=/tmp/mbv-reproducibility_inputs \
           MBV_REPRODUCIBILITY_INPUTS_EPHEMERAL=1 \
    && chmod a+r /opt/mbv/scripts/verify_python_env.py \
    && bash /opt/mbv/scripts/build_portable_env.sh python \
    && python /opt/mbv/scripts/verify_python_env.py \
    && python -m pip check

# Confirm the temporary build inputs were not persisted in the image.
RUN test ! -e /opt/mbv/reproducibility_inputs/conda-linux-64.lock \
    && test ! -e /opt/mbv/reproducibility_inputs/requirements-pip.lock \
    && test ! -e /opt/mbv/reproducibility_inputs/wheels \
    && test ! -e /tmp/mbv-reproducibility_inputs/conda-linux-64.lock \
    && test ! -e /tmp/mbv-reproducibility_inputs/requirements-pip.lock \
    && test ! -e /tmp/mbv-reproducibility_inputs/wheels

WORKDIR /work
CMD ["python", "--version"]
