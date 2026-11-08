#!/usr/bin/env bash
set -o pipefail

REPO=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd
OUTDIR="$REPO/reproducibility_Rcheck"
JOB_ID=${SLURM_JOB_ID:-manual}
TEST_DIR="/tmp/mbv-renv-test-${JOB_ID}"
BOOTSTRAP_LIB="$TEST_DIR/bootstrap-lib"
PROJECT_DIR="$TEST_DIR/project"
TARGET_LIB="$TEST_DIR/renv-library"
CACHE_DIR="$TEST_DIR/cache"

mkdir -p "$PROJECT_DIR" "$BOOTSTRAP_LIB" "$TARGET_LIB" "$CACHE_DIR"
cp "$REPO/renv.lock" "$PROJECT_DIR/renv.lock"
cd "$PROJECT_DIR"
module load conda_R/4.4.x
export R_LIBS_USER="$BOOTSTRAP_LIB"
export RENV_PATHS_CACHE="$CACHE_DIR"
export RENV_CONFIG_AUTO_SNAPSHOT=FALSE

RESTORE_LOG="$OUTDIR/R_restore_6G_job_${JOB_ID}.log"
BOOTSTRAP_LOG="$OUTDIR/R_biocversion_bootstrap_6G_job_${JOB_ID}.log"
POST_LOG="$OUTDIR/R_post_restore_6G_job_${JOB_ID}.log"

{
  echo "Started: $(date -Is)"
  echo "Job: ${JOB_ID}; host: $(hostname); requested memory: 6G"
  echo "R module: conda_R/4.4.x"
  Rscript --vanilla -e 'cat("R version:", R.version.string, "\n"); cat("Bioconductor lock version: 3.20\n")'
  echo "Installing bootstrap renv into isolated bootstrap library"
  Rscript --vanilla -e 'install.packages("renv", repos="https://cran.rstudio.com", lib=Sys.getenv("R_LIBS_USER"))'
} 2>&1 | tee "$RESTORE_LOG"
RENV_BOOTSTRAP_STATUS=${PIPESTATUS[0]}
if [ "$RENV_BOOTSTRAP_STATUS" -ne 0 ]; then
  echo "renv bootstrap failed with exit code $RENV_BOOTSTRAP_STATUS" | tee -a "$RESTORE_LOG"
  exit "$RENV_BOOTSTRAP_STATUS"
fi

# BiocManager's staged load check expects BiocVersion to exist. Seed the exact
# Bioconductor 3.20 package from its pinned repository into the isolated target library.
MBV_TARGET_LIB="$TARGET_LIB" Rscript --vanilla -e 'lib <- Sys.getenv("MBV_TARGET_LIB"); install.packages("BiocVersion", repos=c(BioCsoft="https://bioconductor.org/packages/3.20/bioc"), lib=lib); stopifnot(as.character(packageVersion("BiocVersion", lib.loc=lib)) == "3.20.0"); cat("BiocVersion", as.character(packageVersion("BiocVersion", lib.loc=lib)), "installed in", lib, "\n")' 2>&1 | tee "$BOOTSTRAP_LOG"
BIOCVERSION_STATUS=${PIPESTATUS[0]}

if [ "$BIOCVERSION_STATUS" -eq 0 ]; then
  MBV_TARGET_LIB="$TARGET_LIB" Rscript --vanilla -e 'library(renv, lib.loc=Sys.getenv("R_LIBS_USER")); options(renv.install.timeout=21600); cat("renv install timeout:", getOption("renv.install.timeout"), "seconds\n"); renv::restore(project=getwd(), lockfile="renv.lock", library=Sys.getenv("MBV_TARGET_LIB"), prompt=FALSE, clean=FALSE)' 2>&1 | tee -a "$RESTORE_LOG"
  RESTORE_STATUS=${PIPESTATUS[0]}
else
  RESTORE_STATUS=$BIOCVERSION_STATUS
  echo "BiocVersion bootstrap failed; skipping restore." | tee -a "$RESTORE_LOG"
fi

MBV_TARGET_LIB="$TARGET_LIB" R_LIBS_USER="$BOOTSTRAP_LIB" RENV_PATHS_CACHE="$CACHE_DIR" Rscript --vanilla -e 'library(renv, lib.loc=Sys.getenv("R_LIBS_USER")); lib <- Sys.getenv("MBV_TARGET_LIB"); ip <- installed.packages(lib.loc=lib); cat("R version:", R.version.string, "\n"); cat("Installed in target library:", nrow(ip), "\n"); for (p in c("BiocManager","BiocVersion","Rcpp","Seurat","SingleCellExperiment","spatialLIBD")) cat(p, "installed=", if (p %in% rownames(ip)) ip[p,"Version"] else "NO", "\n", sep=""); status <- renv::status(project=getwd(), library=lib, lockfile="renv.lock", sources=FALSE, cache=FALSE, dev=FALSE); cat("LOCKFILE_SYNCHRONIZED=", isTRUE(status$synchronized), "\n", sep="")' 2>&1 | tee "$POST_LOG"
POST_STATUS=${PIPESTATUS[0]}
echo "Restore exit code: $RESTORE_STATUS; post-check exit code: $POST_STATUS" | tee -a "$POST_LOG"
echo "Finished: $(date -Is)" | tee -a "$RESTORE_LOG"
exit "$RESTORE_STATUS"
