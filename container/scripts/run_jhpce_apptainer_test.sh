#!/usr/bin/env bash
set -euo pipefail

profile=${1:?Usage: run_jhpce_apptainer_test.sh r|python|all}
case "$profile" in
  r|python) profiles=("$profile") ;;
  all) profiles=(r python) ;;
  *) echo "Profile must be r, python, or all" >&2; exit 2 ;;
esac

if [[ -z ${SLURM_JOB_ID:-} ]]; then
  echo "Run this script inside a JHPCE compute allocation." >&2
  echo "Example: srun --mem 10G --pty bash" >&2
  exit 2
fi

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"
module load apptainer/1.5.3
command -v apptainer >/dev/null || {
  echo "Could not load Apptainer 1.5.3" >&2
  exit 127
}

run_id=${MBV_JHPCE_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}
run_root="$repo_root/reproducibility_outputs/apptainer/$run_id"
workdir="/tmp/mbv-apptainer-${SLURM_JOB_ID}"
mkdir -p "$run_root" "$workdir/tmp" "$workdir/cache"
export APPTAINER_TMPDIR="$workdir/tmp"
export APPTAINER_CACHEDIR="$workdir/cache"

cleanup() {
  rc=$?
  if [[ $rc -eq 0 ]]; then
    rm -rf "$workdir"
    printf 'Removed node-local build temp/cache: %s\n' "$workdir" \
      | tee -a "$run_root/cleanup.log"
  else
    printf 'Build temp/cache retained until the Slurm job ends: %s\n' "$workdir" \
      | tee -a "$run_root/cleanup.log"
  fi
}
trap cleanup EXIT

{
  printf 'site=JHPCE\nrun_id_utc=%s\nprofile=%s\n' "$run_id" "$profile"
  printf 'repository=%s\nbranch=%s\ncommit=%s\n' \
    "$repo_root" "$(git branch --show-current)" "$(git rev-parse HEAD)"
  printf 'job_id=%s\nhost=%s\n' "$SLURM_JOB_ID" "$(hostname)"
  printf 'requested_memory=%s\n' "${MBV_JHPCE_MEM_REQUEST:-unknown}"
  scontrol show job "$SLURM_JOB_ID" -o 2>&1 || true
  printf '\nApptainer version: '
  apptainer --version
  printf '\nHost memory:\n'
  free -h
  printf '\nFilesystem capacity:\n'
  df -h "$repo_root" "$workdir"
  printf '\n'
} > "$run_root/session-info.txt"

run_one() {
  local name=$1
  local out="$run_root/${name}-jhpce"
  local definition="$repo_root/container/apptainer/${name}-jhpce.def"
  local image="$out/mbv-$name-jhpce.sif"
  mkdir -p "$out"

  sed "s/^profile=.*/profile=$name/" "$run_root/session-info.txt" > "$out/session-info.txt"
  cat > "$out/commands.sh" <<EOF
#!/usr/bin/env bash
# Run the following allocation request from the JHPCE login node:
srun --mem ${MBV_JHPCE_MEM_REQUEST:-10G} --pty bash
# Then run these commands from the allocated compute node:
cd "$repo_root"
module load apptainer/1.5.3
export MBV_JHPCE_MEM_REQUEST=${MBV_JHPCE_MEM_REQUEST:-10G}
export MBV_JHPCE_RUN_ID=$run_id
export APPTAINER_TMPDIR="$workdir/tmp"
export APPTAINER_CACHEDIR="$workdir/cache"
bash container/scripts/run_jhpce_apptainer_test.sh $name
EOF

  if [[ ! -f $definition ]]; then
    echo "Missing JHPCE definition: $definition" >&2
    return 2
  fi

  cp "$definition" "$out/definition-used.def"
  {
    printf 'profile=%s\nrun_id_utc=%s\nsite=JHPCE\n' "$name" "$run_id"
    printf 'repository=%s\nbranch=%s\ncommit=%s\n' \
      "$repo_root" "$(git branch --show-current)" "$(git rev-parse HEAD)"
    printf 'job_id=%s\nhost=%s\n' "$SLURM_JOB_ID" "$(hostname)"
    printf 'definition=%s\n' "$definition"
    printf 'build_command=apptainer build --fakeroot %q %q\n' "$image" "$definition"
    printf 'test_command=apptainer test %q\n' "$image"
  } > "$out/run-info.txt"
  printf 'requested_memory=%s\n' "${MBV_JHPCE_MEM_REQUEST:-unknown}" >> "$out/run-info.txt"

  local start elapsed rc
  start=$(date +%s)
  if apptainer build --fakeroot "$image" "$definition" > "$out/build.log" 2>&1; then
    elapsed=$(($(date +%s)-start))
    printf 'build_elapsed_seconds=%s\nbuild_status=passed\n' "$elapsed" >> "$out/run-info.txt"
  else
    rc=$?
    elapsed=$(($(date +%s)-start))
    printf 'build_elapsed_seconds=%s\nbuild_status=failed\nbuild_exit_code=%s\n' \
      "$elapsed" "$rc" >> "$out/run-info.txt"
    cat "$out/build.log" >&2
    return "$rc"
  fi

  start=$(date +%s)
  if apptainer test "$image" > "$out/test.log" 2>&1; then
    elapsed=$(($(date +%s)-start))
    printf 'explicit_test_elapsed_seconds=%s\nexplicit_test_status=passed\n' "$elapsed" >> "$out/run-info.txt"
  else
    rc=$?
    elapsed=$(($(date +%s)-start))
    printf 'explicit_test_elapsed_seconds=%s\nexplicit_test_status=failed\nexplicit_test_exit_code=%s\n' \
      "$elapsed" "$rc" >> "$out/run-info.txt"
    cat "$out/test.log" >&2
    return "$rc"
  fi

  sha256sum "$image" > "$out/image.sha256"
  printf 'image_size_bytes=%s\nimage_sha256=%s\n' \
    "$(wc -c < "$image" | tr -d ' ')" \
    "$(awk '{print $1}' "$out/image.sha256")" >> "$out/run-info.txt"
  printf 'JHPCE %s build and test passed. Results: %s\n' "$name" "$out"
  sstat -j "${SLURM_JOB_ID}.${SLURM_STEP_ID:-0}" \
    --format=JobID,MaxRSS,AveRSS,MaxVMSize --parsable2 \
    > "$out/slurm-live-resource-usage.txt" 2>&1 || true
}

overall_status=0
for name in "${profiles[@]}"; do
  if ! run_one "$name"; then
    overall_status=1
  fi
done

exit "$overall_status"
