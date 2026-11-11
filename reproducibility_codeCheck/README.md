# Docker Python AUCell code check

This is a lightweight, script-level check that the `mbv-python:test` Docker image can load the project inputs and run the AUCell analysis code in `docker-aucell.py`. It complements the image smoke test, which checks the installed Python environment but does not execute this project analysis script.

## Why the loom input is subset to 10%

The full loom input contains data for 119 samples and caused the AUCell process to be killed when run with Docker Desktop's 8 GiB memory limit (`OOMKilled=true`, exit status 137). The script therefore writes a separate loom file containing the first 10% of cells, then runs AUCell on that smaller file. It uses one worker as an additional memory precaution.

The subset is intended to check that the code runs end to end on a smaller input. It is not a full-data pipeline run, a validation of full-data memory requirements, or a substitute for checking scientific results on all cells. The subset file is created alongside the source loom as `*_first10pct.loom`; the original loom is left unchanged. Loom files are compressed, so the subset's disk size may not be exactly 10% of the original.

## Requirements and inputs

- Run the command from the cloned MBv project root.
- Docker Desktop must be running, and the `mbv-python:test` image must already be built locally for `linux/amd64`.
- The project test inputs must be present:
  - `reproducibility_codeCheck/data/expr_loom/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts.loom`
  - `reproducibility_codeCheck/data/spe-n119_13162-no-lowUMI_modules_DEG-subset-refined.dat`
- The `test_results/docker/` directory must exist for the log and AUCell output.

## Run the check

In zsh, from the project root:

```zsh
# Add Docker's CLI to PATH if this shell cannot find the `docker` command.
export PATH="/Applications/Docker.app/Contents/Resources/bin:$PATH"
docker info

# Preserve the Python output in a log and report the container's exit status.
set -o pipefail
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" \
  -w /work \
  mbv-python:test \
  python3 -u reproducibility_codeCheck/docker-aucell.py \
  2>&1 | tee test_results/docker/aucell.log
print "exit status: $?"
```

The bind mount makes the project files available at `/work` in the container and keeps generated files on the host checkout. `tee` replaces `test_results/docker/aucell.log` on each run; use `tee -a` if you intentionally want to append instead.

## Outputs and pass/fail

- `test_results/docker/aucell.log` captures standard output and errors from the run.
- `test_results/docker/spe-n119_seurat-label-Micro.Vasc_13162-genes_no-lowUMI_logcounts_first10pct_modules-DEG-subset-refined_AUCell-fixed-thold-05.loom` is the AUCell output loom.
- The generated input subset is written next to the original loom with the suffix `_first10pct.loom`.
- Exit status `0` indicates the script completed without a reported error. A nonzero status means the run failed; inspect the log for the traceback or termination details. In particular, status `137` together with Docker reporting `OOMKilled=true` means the process exceeded its memory limit.

This check exercises only the Docker Python image and the reduced AUCell test input. It does not test the R image, Apptainer images, or a full-data run.
