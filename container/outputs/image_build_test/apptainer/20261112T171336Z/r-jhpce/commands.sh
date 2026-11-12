#!/usr/bin/env bash
# Run the following allocation request from the JHPCE login node:
srun --mem 10G --pty bash
# Then run these commands from the allocated compute node:
cd "/users/jthompso/spatialDLPFC_mdd_bpd"
module load apptainer/1.5.3
export MBV_JHPCE_MEM_REQUEST=10G
export MBV_JHPCE_RUN_ID=20261112T171336Z
export APPTAINER_TMPDIR="/tmp/mbv-apptainer-37038242/tmp"
export APPTAINER_CACHEDIR="/tmp/mbv-apptainer-37038242/cache"
bash container/02_image_build_test/scripts/run_apptainer_jhpce_build_smoke_test.sh r
