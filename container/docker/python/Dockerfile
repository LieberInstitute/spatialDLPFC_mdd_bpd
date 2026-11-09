FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="MBv Python environment" \
      org.opencontainers.image.version="0.1-test" \
      org.opencontainers.image.description="pySCENIC Conda environment for the MBv pipeline"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential bash ca-certificates curl bzip2 procps \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/mbv/reproducibility_inputs /opt/mbv/scripts

COPY reproducibility_inputs/conda-linux-64.lock /opt/mbv/reproducibility_inputs/conda-linux-64.lock
COPY reproducibility_inputs/requirements-pip.lock /opt/mbv/reproducibility_inputs/requirements-pip.lock
COPY reproducibility_inputs/wheels/ /opt/mbv/reproducibility_inputs/wheels/
COPY container/scripts/verify_python_env.py /opt/mbv/scripts/verify_python_env.py
COPY container/scripts/build_portable_env.sh /opt/mbv/scripts/build_portable_env.sh

ENV PATH=/opt/conda/envs/pyscenic/bin:/opt/conda/bin:$PATH

RUN chmod a+r /opt/mbv/scripts/verify_python_env.py \
    && bash /opt/mbv/scripts/build_portable_env.sh python \
    && python /opt/mbv/scripts/verify_python_env.py \
    && python -m pip check

WORKDIR /work
CMD ["python", "--version"]
