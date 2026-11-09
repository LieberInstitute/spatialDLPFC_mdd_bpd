#!/usr/bin/env bash
set -eo pipefail

REPO=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd
OUTDIR="$REPO/reproducibility_apptainer"
JOB_ID=${SLURM_JOB_ID:-manual}
WORKDIR="/tmp/mbv-apptainer-test-${JOB_ID}"
SIF="$WORKDIR/MBv-pipeline.sif"

mkdir -p "$OUTDIR" "$WORKDIR/tmp" "$WORKDIR/cache"
export APPTAINER_TMPDIR="$WORKDIR/tmp"
export APPTAINER_CACHEDIR="$WORKDIR/cache"

module load apptainer/1.5.3

cat > "$OUTDIR/apptainer_test_commands_${JOB_ID}.txt" <<EOF
JHPCE revision branch, commit: $(git -C "$REPO" rev-parse --short HEAD)
Compute node request: srun --mem 10G --pty bash
Apptainer module: apptainer/1.5.3
Build: apptainer build --fakeroot $SIF $REPO/container/MBv-pipeline.def
Container checks: apptainer test $SIF
Post-build smoke: Rscript --version; python --version; python -m pip check
Temporary SIF and build cache: $WORKDIR (removed after successful checks)
EOF

apptainer --version | tee "$OUTDIR/apptainer_version_${JOB_ID}.txt"
apptainer build --fakeroot "$SIF" "$REPO/container/MBv-pipeline.def" \
    2>&1 | tee "$OUTDIR/apptainer_build_${JOB_ID}.log"
apptainer test "$SIF" 2>&1 | tee "$OUTDIR/apptainer_test_${JOB_ID}.log"
apptainer exec "$SIF" Rscript --version \
    2>&1 | tee "$OUTDIR/apptainer_R_version_${JOB_ID}.txt"
apptainer exec "$SIF" python --version \
    2>&1 | tee "$OUTDIR/apptainer_python_version_${JOB_ID}.txt"
apptainer exec "$SIF" python -m pip check \
    2>&1 | tee "$OUTDIR/apptainer_pip_check_${JOB_ID}.txt"
sha256sum "$SIF" | tee "$OUTDIR/apptainer_sif_sha256_${JOB_ID}.txt"

# Keep test records on the project filesystem; discard only temporary build data.
rm -rf "$WORKDIR"
echo "Temporary build directory removed: $WORKDIR" \
    | tee "$OUTDIR/apptainer_cleanup_${JOB_ID}.txt"
